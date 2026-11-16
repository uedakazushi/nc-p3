"""Eight-parameter product-net chart for four quadrics through a fixed point."""
from itertools import combinations_with_replacement,product
import sympy as s
from frobenius_tensor import digits
import generic_boundary_field as g

NAMES=('al','be','ga','u','v','a','b','c')


def terms(K=None):
    variables=s.symbols(NAMES)
    al,be,ga,u,v,a,b,c=(variables if K is None else [K.from_sympy(x) for x in variables])
    one=s.S.One if K is None else K.one
    return [{(0,0):one,(0,3):-one,(0,1):-al},
            {(1,1):one,(1,3):-one,(1,2):-be},
            {(2,2):one,(2,3):-one,(0,2):-ga},
            {(0,3):one,(1,3):u,(2,3):v,(0,1):a,(0,2):b,(1,2):c}]


def tensor(K,polynomial=False):
    """Solve the 20-coordinate cubic frame; optionally clear return-basis denominators."""
    monomials=list(combinations_with_replacement(range(4),3));lookup={m:i for i,m in enumerate(monomials)}
    rows=[]
    for last in range(4):
        for quadric in terms(K):
            row={}
            for pair,value in quadric.items():
                key=lookup[tuple(sorted(pair+(last,)))]
                row[key]=row.get(key,K.zero)+value
            rows.append({i:x for i,x in row.items() if x})
    free=[(0,1,2),(1,3,3),(2,3,3),(3,3,3)]
    original=[dict(row) for row in rows]
    remaining=rows+[{lookup[m]:K.one,20+j:K.one} for j,m in enumerate(free)]
    pivots=[]
    while remaining:
        selected=None;best=None
        for r,row in enumerate(remaining):
            for column,value in row.items():
                if column>=20:continue
                score=(len(value.n)+len(value.d),value.n.total_degree()+value.d.total_degree(),len(row))
                if best is None or score<best:best,selected=score,(r,column)
        assert selected is not None
        r,column=selected;row=remaining.pop(r);scale=row[column]
        row={i:x/scale for i,x in row.items()}
        updated=[]
        for other in remaining:
            factor=other.get(column,K.zero)
            if factor:
                other=dict(other)
                for i,value in row.items():
                    total=other.get(i,K.zero)-factor*value
                    if total:other[i]=total
                    else:other.pop(i,None)
            if other:updated.append(other)
        remaining=updated;pivots.append((column,row))
    assert len(pivots)==20
    basis={}
    for column,row in reversed(pivots):
        basis[column]=[row.get(20+j,K.zero)-sum((value*basis[i][j]
            for i,value in row.items() if i<20 and i!=column),K.zero) for j in range(4)]
    for row in original:
        for j in range(4):assert sum((x*basis[i][j] for i,x in row.items()),K.zero)==K.zero
    for i,m in enumerate(free):assert basis[lookup[m]]==[K.one if i==j else K.zero for j in range(4)]
    denominators=[]
    for j in range(4):
        denominator=K.one.d
        if polynomial:
            for vector in basis.values():
                other=vector[j].d
                denominator=(denominator*other)//denominator.gcd(other)
        denominators.append(type(K.one)(denominator))
    result={}
    for encoded in range(256):
        a,b,j,c=digits(encoded)
        value=basis[lookup[tuple(sorted((a,b,c)))]] [j]*denominators[j]
        if value:result[encoded]=value
    if polynomial:assert all(x.d==1 for x in result.values())
    return result


def relation_windows(w,K,normalize_r0=False):
    al,be,ga,u,v,a,b,c=[K.from_sympy(x) for x in s.symbols(NAMES)]
    raw=[]
    for free in range(6):
        t=[K.one if i==free else K.zero for i in range(6)]
        r03=-(u*t[3]+v*t[4]+a*t[0]+b*t[1]+c*t[2])
        row={(0,0):r03+al*t[0],(1,1):t[3]+be*t[2],(2,2):t[4]+ga*t[1],(3,3):t[5]}
        for pair,value in [((0,1),t[0]),((0,2),t[1]),((1,2),t[2]),((1,3),t[3]),((2,3),t[4]),((0,3),r03)]:
            row[pair]=value;row[pair[::-1]]=value
        raw.append({4*i+j:x for (i,j),x in row.items() if x})
    if normalize_r0:
        normalized=g.echelon(raw);assert len(normalized)==6
        r0=[{divmod(i,4):x for i,x in row.items()} for row in normalized.values()]
    else:r0=[{divmod(i,4):x for i,x in row.items()} for row in raw]
    r1=[]
    for index,(first,last) in enumerate([(0,1),(0,2),(1,2),(1,3),(2,3),(3,3)]):
        row={}
        for b,j in product(range(4),repeat=2):
            value=w.get(64*first+16*b+4*j+last,K.zero)
            if value:row[b,j]=value
        unit_coordinate=[(2,0),(1,0),(0,0),(3,1),(3,2),(3,3)][index]
        scale=row[unit_coordinate]
        r1.append({pair:x/scale for pair,x in row.items()})
    r2=[{(b,a):x for (a,b),x in row.items()} for row in r1]
    return [r0,r1,r2,r0]
