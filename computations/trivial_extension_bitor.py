"""One-return two-sided Koszul complex for the four-layer quadric truncation.

The caller supplies a cubic frame, a six-vector B2 dual frame, and four
commutative quadrics. The complex computes the internal-degree-five,
vertex-(position1,position2) summand of Tor^{A^e}(S^e, DA[degree4]).
"""
from itertools import combinations


def relation_frame(quadrics,K):
    result=[]
    for a,b in combinations(range(4),2):
        result.append({(a,b):K.one,(b,a):-K.one})
    result.extend(quadrics)
    assert len(result)==10
    return result


def one_return_complex(w,b2_dual,quadrics,K):
    """Return sparse differentials d_h:C_h -> C_(h-1), h=2,3,4."""
    zero,one=K.zero,K.one
    q2=relation_frame(quadrics,K)
    qdims=(1,4,10);bdims=(1,4,6,4)
    free2=[(0,1),(0,2),(1,2),(1,3),(2,3),(3,3)]
    # C_i is dual to multiplication by e_i in the commutative ring B.
    contractions={}
    for r in (1,2,3):
        for i in range(4):
            for j in range(bdims[r]):
                if r==1:column={0:one} if i==j else {}
                elif r==2:column={a:b2_dual[j].get((i,a),zero) for a in range(4)}
                else:column={a:w.get(64*i+16*b+4*j+c,zero) for a,(b,c) in enumerate(free2)}
                contractions[r,i,j]={a:x for a,x in column.items() if x}
    bases={}
    for h in range(1,5):
        bases[h]=[(p,h-p,a,j,b) for p in range(3) if 0<=h-p<=2
                  for a in range(qdims[p]) for j in range(bdims[h-1])
                  for b in range(qdims[h-p])]
    indexes={h:{v:i for i,v in enumerate(base)} for h,base in bases.items()}
    differentials={}
    for h in (2,3,4):
        columns=[];r=h-1
        for p,q,a,j,b in bases[h]:
            column={}
            def add(key,value):
                row=indexes[h-1][key];total=column.get(row,zero)+value
                if total:column[row]=total
                else:column.pop(row,None)
            if p:
                pieces={(0,a):one} if p==1 else q2[a]
                for (prefix,arrow),coefficient in pieces.items():
                    for m,value in contractions[r,arrow,j].items():
                        add((p-1,q,prefix,m,b),coefficient*value)
            if q:
                pieces={(b,0):one} if q==1 else q2[b]
                sign=one if p%2==0 else -one
                for (arrow,suffix),coefficient in pieces.items():
                    for m,value in contractions[r,arrow,j].items():
                        add((p,q-1,a,m,suffix),sign*coefficient*value)
            columns.append(column)
        differentials[h]=columns
    assert [len(bases[h]) for h in range(1,5)]==[8,144,480,400]
    return differentials,bases
