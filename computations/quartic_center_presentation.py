"""Certify a Hilbert series by independent relations over five quartic variables."""
from collections import defaultdict
from itertools import product

from verify_klein_center_lift import add_multiple

K = range(4)
# Coefficients of (1-t^4)^5/(1-t)^4.
NUMERATOR = (1,4,10,20,30,36,34,20,0,-20,-34,-36,-30,-20,-10,-4,-1)


def lifted_rhs(rows,rhs,pivots,independent):
    """Apply the selected square minor to the lower-word contributions."""
    if not any(rhs):
        return {j:{} for j in pivots}
    triangular = {}
    values = {}
    for j,index in zip(pivots,independent,strict=True):
        row = rows[index].copy()
        value = rhs[index].copy()
        for col,prior in triangular.items():
            if col in row:
                factor = -row[col]
                add_multiple(row,prior,factor)
                add_multiple(value,values[col],factor)
        divisor = row[j]
        triangular[j] = {k:v/divisor for k,v in row.items()}
        values[j] = {k:v/divisor for k,v in value.items()}
    for j in reversed(triangular):
        for col,c in triangular[j].items():
            if col!=j and col in triangular:
                add_multiple(values[j],values[col],-c)
    return values


def certify_presentation(windows,sections,start,constants,one,pivot,relation_basis,profile):
    """Lift the full-center quotient basis and evaluate homogeneous syzygies at constants."""
    assert any(constants) and len(constants)==5 and profile[-1]==0
    words_by_color = {(0,h):[()] if h==0 else [] for h in K}
    arrows = {}

    def append(vector,i):
        out = {}
        for word,c in vector.items():
            add_multiple(out,arrows[word,i],c)
        return out

    def apply(vector,word):
        for i in word:
            vector = append(vector,i)
        return vector

    for n in range(1,len(profile)):
        for h in K:
            words = [word+(i,) for i in K for word in words_by_color[n-1,h^i]]
            index = {w:j for j,w in enumerate(words)}
            rows = []
            rhs = []
            def equation(terms,target):
                row = {}
                lower = target.copy()
                for vector,i,c in terms:
                    for word,v in vector.items():
                        if len(word)!=n-1:
                            add_multiple(lower,arrows[word,i],-c*v)
                            continue
                        j = index[word+(i,)]
                        val = row.get(j,0)+c*v
                        if val:
                            row[j] = val
                        else:
                            row.pop(j,None)
                rows.append(row)
                rhs.append(lower)
            if n>=2:
                for g,rule in windows[(start+n-2)%4]:
                    for word in words_by_color[n-2,h^g]:
                        equation(((arrows[word,i],j,c) for (i,j),c in rule.items()),{})
            if n==4:
                for k,(g,section) in enumerate(sections[start]):
                    if g==h:
                        equation(((apply({():one},w[:-1]),w[-1],c) for w,c in section.items()),
                                 {():constants[k]} if constants[k] else {})
            pivots,independent = pivot(rows)
            lower = lifted_rhs(rows,rhs,pivots,independent)
            words_by_color[n,h] = [w for j,w in enumerate(words) if j not in pivots]
            for j,w in enumerate(words):
                if j in pivots:
                    arrows[w[:-1],w[-1]] = lower[j]|{words[k]:-v for k,v in pivots[j].items() if k!=j}
                else:
                    arrows[w[:-1],w[-1]] = {w:one}
        assert sum(len(words_by_color[n,h]) for h in K)==profile[n]

    initial = defaultdict(list)
    for (d,_),words in words_by_color.items():
        for word in words:
            for _,rule in windows[(start+d)%4]:
                row = {}
                for pair,c in rule.items():
                    add_multiple(row,apply({word:one},pair),c)
                if row:
                    initial[d+2].append(row)
    for k,(_,section) in enumerate(sections[start]):
        row = {():-constants[k]} if constants[k] else {}
        for word,c in section.items():
            add_multiple(row,apply({():one},word),c)
        if row:
            initial[4].append(row)

    beta = [(profile[n] if n<len(profile) else 0)-v for n,v in enumerate(NUMERATOR)]
    assert len(profile)<=len(beta) and all(n>=0 for n in beta)
    spans = [[]]
    for n in range(1,len(beta)):
        # Multiplication by any central variable with nonzero evaluation
        # identifies every previously certified relation four degrees lower.
        rows = list(initial[n])
        if n>=4:
            rows.extend(spans[n-4])
        rows.extend(append(row,i) for row,i in product(spans[-1],K))
        spanning = relation_basis(rows)
        assert len(spanning)>=sum(beta[n::-4]),(start,n,len(spanning),sum(beta[n::-4]))
        spans.append(list(spanning.values()))
