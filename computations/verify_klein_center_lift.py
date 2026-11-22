"""Exact central quartics over the two-variable slice of the second completion."""
from functools import cache, reduce
from itertools import product
from operator import xor
import sympy as sp
from flint_rational_function import rational_function_field
from verify_klein_period_four_center import K, artinian_quotient, commuting_space, paths, relations_from_gram
from verify_klein_two_doublets import gram


def parameters(b,c,v,w):
    z = (1-c)/(b-1)
    ratio = (1-b*c)/(v*w)
    s = [[1+v,b],[c,1+w]]
    second = [[1+ratio*v,b],[c,1+ratio*w]]
    return z,s,second


def transpose(mat):
    return list(zip(*mat))


def right_inverse_product(m,n):
    (a,b),(c,d)=n
    det=a*d-b*c
    return [[(x*d-y*c)/det,(y*a-x*b)/det] for x,y in m]


def field_relations(b,c,v,w):
    z,s,second=parameters(b,c,v,w)
    x,y=[[1,1],[1,z]],[[1,1],[z,1]]
    graphs=[(right_inverse_product(y,s),right_inverse_product(transpose(x),transpose(second))),
            (right_inverse_product(transpose(x),transpose(s)),right_inverse_product(y,second))]
    relations=[]
    for side,thick in enumerate(graphs):
        rows=[]
        for h,graph in enumerate(thick):
            for i in range(2):
                coeff=[int(i==j) for j in range(2)]+[graph[j][i] for j in range(2)]
                rows.append((h,{(j,j^h):a for j,a in enumerate(coeff) if a}))
        for h in (2,3):
            rows.append((h,{(i,i^h):z if i==3 and (i^h^side)==1 else 1 for i in K}))
        relations.extend([rows,[(h,{(j,i):a for (i,j),a in row.items()}) for h,row in rows]])
    return relations


def add_multiple(target,row,factor):
    if factor:
        for j,a in row.items():
            value=target.get(j,0)+factor*a
            if value:
                target[j]=value
            else:
                target.pop(j,None)


def echelon(rows):
    pivots={}
    for original in rows:
        row=original.copy()
        for j in sorted(pivots):
            if j in row:
                add_multiple(row,pivots[j],-row[j])
        if row:
            j=min(row)
            pivots[j]={k:v/row[j] for k,v in row.items()}
    for j in sorted(pivots,reverse=True):
        for k in list(pivots[j]):
            if k>j and k in pivots:
                add_multiple(pivots[j],pivots[k],-pivots[j][k])
    return pivots


def field_paths(relations,one,reducer=echelon):
    basis={(p,0,h):[()] if h==0 else [] for p,h in product(K,repeat=2)}
    prepend={}
    @cache
    def image(p,word):
        if not word:
            return {0:one}
        h=reduce(xor,word)
        projection=prepend[p,len(word),h,word[0]]
        out={}
        for j,v in image((p+1)%4,word[1:]).items():
            add_multiple(out,projection[j],v)
        return out
    for n in range(1,6):
        for p,h in product(K,repeat=2):
            words=[(i,)+w for i in K for w in basis[(p+1)%4,n-1,h^i]]
            off=[sum(w[0]<i for w in words) for i in range(5)]
            rows=[]
            if n>=2:
                for g,rel in relations[p]:
                    for tail in basis[(p+2)%4,n-2,h^g]:
                        row={}
                        for (i,j),a in rel.items():
                            for k,v in image((p+1)%4,(j,)+tail).items():
                                col=off[i]+k
                                row[col]=row.get(col,0)+a*v
                        rows.append({k:v for k,v in row.items() if v})
            piv=reducer(rows)
            free=[j for j in range(len(words)) if j not in piv]
            ids={j:i for i,j in enumerate(free)}
            basis[p,n,h]=[words[i] for i in free]
            projection=[{ids[j]:one} if j in ids else
                        {ids[k]:-v for k,v in piv[j].items() if k!=j} for j in range(len(words))]
            for i in K:
                prepend[p,n,h,i]=projection[off[i]:off[i+1]]
    return basis,image


def central_equations(basis,image,h):
    dims=[len(basis[p,4,h]) for p in K]
    offsets=[sum(dims[:p]) for p in range(5)]
    rows=[]
    for p,i in product(K,repeat=2):
        nxt=(p+1)%4
        block=[{} for _ in basis[p,5,h^i]]
        for j,w in enumerate(basis[p,4,h]):
            for k,v in image(p,w+(i,)).items():
                block[k][offsets[p]+j]=v
        for j,w in enumerate(basis[nxt,4,h]):
            for k,v in image(p,(i,)+w).items():
                block[k][offsets[nxt]+j]=-v
        rows.extend(block)
    return rows,offsets[-1]


def sample():
    values=tuple(map(sp.Integer,(2,3,4,6)))
    z,s,t=parameters(*values)
    w=gram(sp.ones(2),sp.Matrix([[1,1],[1,z]]),sp.Matrix(s),sp.Matrix(t))
    rel=relations_from_gram(w)
    for original,sparse in zip(rel,field_relations(*values)):
        a=sp.Matrix([r for _,r in original])
        b=sp.Matrix([[r.get(ij,0) for ij in product(K,repeat=2)] for _,r in sparse])
        assert a.col_join(b).rank()==b.rank()==6
    basis,image=paths(rel,5)
    assert [commuting_space(basis,image,h,0).cols for h in K]==[3,2,0,0]
    artinian_quotient(rel,basis,image)


if __name__=='__main__':
    sample()
    field=rational_function_field(('b','c'))
    b,c=(field.from_sympy(t) for t in sp.symbols('b c'))
    basis,image=field_paths(field_relations(b,c,4*field.one,6*field.one),field.one)
    for h,expected in [(0,3),(1,2)]:
        rows,size=central_equations(basis,image,h)
        assert size-len(echelon(rows))==expected
        print('Rational function field: character',h,'quartic center dimension',expected,flush=True)
