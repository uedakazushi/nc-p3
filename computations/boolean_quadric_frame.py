"""Cubic frames for quadrics through exactly one Boolean intersection point."""
from itertools import product
import sympy as s
from frobenius_tensor import digits


def data(field):
    u,v,a,b,c=(field.from_sympy(x) for x in s.symbols('u v a b c'))
    one=field.one
    values={(1,0,0):one,(0,1,0):u,(0,0,1):v,
            (1,1,0):a,(1,0,1):b,(0,1,1):c,
            (1,1,1):a+b+c-one-u-v}
    quadrics={(0,3):one,(1,3):u,(2,3):v,
              (0,1):a-one-u,(0,2):b-one-v,(1,2):c-u-v}
    return (u,v,a,b,c),values,quadrics


def evaluation_tensor(field):
    _,values,_=data(field)
    result={}
    for encoded in range(256):
        a,b,j,c=digits(encoded)
        if j==3:
            value=field.one if (a,b,c)==(3,3,3) else field.zero
        else:
            value=field.zero
            for point,q in values.items():
                if point[j]!=1:continue
                if all(i==3 or point[i] for i in (a,b,c)):
                    sign=(-1)**(sum(point)-1)
                    value+=field.from_sympy(s.Integer(sign))/q
        if value:result[encoded]=value
    return result


def normalized_tensor(field):
    (u,v,a,b,c),_,quadrics=data(field)
    z,o=field.zero,field.one
    # Unknown moments m0,m01,m02,m12; free moments m012,m1,m2,m_empty.
    matrix=[[o,a-o,b-o,z],
            [z,a-u,z,c-u],
            [z,z,b-v,c-v],
            [o,a-o-u,b-o-v,c-u-v]]
    rhs=[[-(c-u-v),z,z,z],[-(b-o-v),-u,z,z],
         [-(a-o-u),z,-v,z],[z,-u,-v,z]]
    augmented=[row+tail for row,tail in zip(matrix,rhs)]
    for pivot in range(4):
        chosen=next(i for i in range(pivot,4) if augmented[i][pivot])
        augmented[pivot],augmented[chosen]=augmented[chosen],augmented[pivot]
        scale=augmented[pivot][pivot]
        augmented[pivot]=[x/scale for x in augmented[pivot]]
        for i in range(4):
            if i==pivot:continue
            scale=augmented[i][pivot]
            if scale:augmented[i]=[x-scale*y for x,y in zip(augmented[i],augmented[pivot])]
    moments={frozenset((0,)):augmented[0][4:],frozenset((0,1)):augmented[1][4:],
             frozenset((0,2)):augmented[2][4:],frozenset((1,2)):augmented[3][4:],
             frozenset((0,1,2)):[o,z,z,z],frozenset((1,)):[z,o,z,z],
             frozenset((2,)):[z,z,o,z],frozenset():[z,z,z,o]}
    def coefficient(word,j):return moments[frozenset(i for i in word if i!=3)][j]
    for j,a in product(range(4),repeat=2):
        assert sum((q*coefficient(pair+(a,),j) for pair,q in quadrics.items()),z)==z
    result={}
    for encoded in range(256):
        a,b,j,c=digits(encoded);value=coefficient((a,b,c),j)
        if value:result[encoded]=value
    return result


def relation_windows(w,field,clear_normalization=True):
    """Generic relation bases without an expensive dense flattening RREF.

The middle windows are contractions by a fixed six-element basis of B2.
Their generic independence follows already from the checked u=v=1,
(a,b,c)=(2,2,2) point. No claim of this basis at every special point is made.
"""
    import generic_boundary_field as g
    (u,v,a,b,c),_,quadrics=data(field)
    d01,d02,d12=quadrics[0,1],quadrics[0,2],quadrics[1,2]
    def si(i):return {(i,i):field.one,(i,3):field.one,(3,i):field.one}
    s0=si(0)
    raw=[]
    for i,j,coefficient in [(0,1,d01),(0,2,d02),(1,2,d12)]:
        relation={(i,j):field.one,(j,i):field.one}
        for pair,value in s0.items():relation[pair]=-coefficient*value
        raw.append({p:x for p,x in relation.items() if x})
    for i,coefficient in [(1,u),(2,v)]:
        relation=si(i)
        for pair,value in s0.items():relation[pair]=-coefficient*value
        raw.append(relation)
    raw.append({(3,3):field.one})
    normalized=g.echelon([{4*i+j:x for (i,j),x in row.items()} for row in raw])
    assert len(normalized)==6
    r0=[{divmod(j,4):x for j,x in row.items()} for row in normalized.values()]
    determinant=(-(a-u)*(b-v)*(c-u-v)-v*(a-u)*(c-v)-u*(b-v)*(c-u)
                 if clear_normalization else field.one)
    r1=[]
    for first,last in [(0,1),(0,2),(1,2),(1,3),(2,3),(3,3)]:
        row={}
        for b,j in product(range(4),repeat=2):
            value=w.get(64*first+16*b+4*j+last,field.zero)*determinant
            if value:row[b,j]=value
        r1.append(row)
    r2=[{(b,a):x for (a,b),x in row.items()} for row in r1]
    return [r0,r1,r2,r0]


def polynomial_tensor(field):
    """Rescale the three face vectors by units to avoid all denominators."""
    _,values,_=data(field)
    scales=[]
    for i in range(3):
        scale=field.one
        for point,q in values.items():
            if point[i]:scale*=q
        scales.append(scale)
    scales.append(field.one)
    result={i:value*scales[digits(i)[2]] for i,value in evaluation_tensor(field).items()}
    assert all(value.d==1 for value in result.values())
    return result
