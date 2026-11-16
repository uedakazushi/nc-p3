"""Cubic frame after eliminating one h from x_i*h=g_i(x).

Only a ten-coordinate ternary cubic system is solved. No homology rank
is inferred from this reduction.
"""
from itertools import combinations,combinations_with_replacement
import generic_boundary_field as field
from product_net_boundary import terms
from frobenius_tensor import digits


def ternary_data(K):
    quadrics=terms(K)
    g=[]
    for i in range(3):
        g.append({pair:value for pair,value in quadrics[i].items() if 3 not in pair and value})
    q={pair:value for pair,value in quadrics[3].items() if value}
    coefficients=[q.get((i,3),K.zero) for i in range(3)]
    for i,c in enumerate(coefficients):
        for pair,value in quadrics[i].items():
            total=q.get(pair,K.zero)+c*value
            if total:q[pair]=total
            else:q.pop(pair,None)
    assert not any(3 in pair for pair in q)
    return g,q


def _times_x(quadric,i):
    return {tuple(sorted(pair+(i,))):value for pair,value in quadric.items()}


def equations(K):
    g,q=ternary_data(K);rows=[]
    for i in range(3):rows.append(_times_x(q,i))
    for i,j in combinations(range(3),2):
        row=_times_x(g[i],j)
        for word,value in _times_x(g[j],i).items():
            total=row.get(word,K.zero)-value
            if total:row[word]=total
            else:row.pop(word,None)
        rows.append(row)
    row={}
    for (a,b),coefficient in q.items():
        for word,value in _times_x(g[a],b).items():
            total=row.get(word,K.zero)+coefficient*value
            if total:row[word]=total
            else:row.pop(word,None)
    rows.append(row)
    assert len(rows)==7
    return rows


def tensor(K):
    g,q=ternary_data(K)
    mons=list(combinations_with_replacement(range(3),3));lookup={m:i for i,m in enumerate(mons)}
    rows=[{lookup[m]:x for m,x in row.items()} for row in equations(K)]
    normal=[(0,1,2),(1,1,1),(2,2,2)]
    rows.extend({lookup[m]:K.one,10+j:K.one} for j,m in enumerate(normal))
    # RHS variables are kept out of the pivot selection.
    remaining=rows;history=[]
    while remaining:
        candidates=[((len(value.n)+len(value.d),value.n.total_degree()+value.d.total_degree(),len(row)),r,i)
                    for r,row in enumerate(remaining) for i,value in row.items() if i<10 and value]
        assert candidates,'The selected ternary normalization chart is singular.'
        _,r,i=min(candidates);row=remaining.pop(r);scale=row[i]
        row={j:x/scale for j,x in row.items()}
        other_rows=[]
        for other in remaining:
            factor=other.get(i,K.zero)
            if factor:
                other=dict(other)
                for j,value in row.items():
                    total=other.get(j,K.zero)-factor*value
                    if total:other[j]=total
                    else:other.pop(j,None)
            if other:other_rows.append(other)
        remaining=other_rows;history.append((i,row))
    assert len(history)==10
    basis={}
    for i,row in reversed(history):
        basis[i]=[row.get(10+j,K.zero)-sum((value*basis[c][j]
            for c,value in row.items() if c<10 and c!=i),K.zero) for j in range(3)]
    for row in equations(K):
        for j in range(3):assert not sum((x*basis[lookup[m]][j] for m,x in row.items()),K.zero)
    values={}
    def value(word,j):
        word=tuple(sorted(word))
        if j==3:return K.one if word==(3,3,3) else K.zero
        if word in values:return values[word][j]
        if 3 not in word:answer=basis[lookup[word]]
        elif word==(3,3,3):answer=[K.zero]*3
        else:
            a=next(i for i in word if i!=3);tail=list(word);tail.remove(a);tail.remove(3)
            answer=[sum((coefficient*value(pair+tuple(tail),r) for pair,coefficient in g[a].items()),K.zero) for r in range(3)]
        values[word]=answer
        return answer[j]
    result={}
    for n in range(256):
        a,b,j,c=digits(n);x=value((a,b,c),j)
        if x:result[n]=x
    # Check the original four-variable cubic annihilation equations.
    for quadric in terms(K):
        for last in range(4):
            for j in range(4):
                assert not sum((coefficient*value(pair+(last,),j) for pair,coefficient in quadric.items()),K.zero)
    return result


def adapted_relations(K):
    """Ten relations with separate left/right h-leading terms."""
    g,q=ternary_data(K);rows=[]
    for i in range(3):
        for pair in ((i,3),(3,i)):
            row={pair:K.one}
            for word,x in g[i].items():row[word]=-x
            rows.append(row)
    for i,j in combinations(range(3),2):rows.append({(i,j):K.one,(j,i):-K.one})
    rows.append(q)
    return rows
