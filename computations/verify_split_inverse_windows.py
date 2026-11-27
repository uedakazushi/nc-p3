"""Complementary inverses for the two alternating sectors of a two-plane splitting."""
from itertools import product

import sympy as sp

BITS = tuple(product(range(2),repeat=4))
PAIRS = ((3,2,1,0), (1,0,3,2))
PATTERNS = ((0,0,0,0),(0,0,1,1),(0,1,0,1),(0,1,1,0))


def positions(i,j,start):
    q = (i//2,i%2,j//2,j%2)
    return q[-start:]+q[:-start] if start else q


def cut(tensor,start):
    return sp.Matrix(4,4,lambda i,j: tensor[positions(i,j,start)])


def reciprocal(tensor,start):
    matrix = cut(tensor,start).inv().T
    return {positions(i,j,start):matrix[i,j] for i,j in product(range(4),repeat=2)}


def transport(tensor, maps, pairing):
    return {q:sum(v*sp.prod(maps[i][q[pairing[i]],j[i]] for i in range(4))
                  for j,v in tensor.items()) for q in BITS}


def borders(maps,pairing,start,q):
    # Row side of the cut uses alternating spaces, column side uses their complementary spaces.
    first = {start,(start+1)%4}
    return (sp.prod(maps[i][q[pairing[i]],q[i]] for i in first),
            sp.prod(maps[i][q[pairing[i]],q[i]] for i in range(4) if i not in first))


def schur_transport():
    p = {q:sp.Integer(3+sum((2*i+1)*b for i,b in enumerate(q))**3+int(q[0]==q[1]==q[2]==q[3]))
         for q in BITS}
    maps = [sp.Matrix([[i+1,2],[1,3+i]]) for i in range(4)]
    for start,pairing in enumerate(PAIRS):
        q = transport(reciprocal(p,start), maps, pairing)
        boundary = {w:borders(maps,pairing,start,w) for w in BITS}
        u = cut({w:v[0] for w,v in boundary.items()},start)
        v = cut({w:v[1] for w,v in boundary.items()},start)
        assert v*cut(p,start).inv()*u == cut(q,start)


def bell(x):
    a,b,c,d = x
    values = ((a+b)/2,(a-b)/2,(c+d)/2,(c-d)/2)
    result = dict.fromkeys(BITS,0)
    for q,value in zip(PATTERNS,values):
        result[q] = result[tuple(1-i for i in q)] = value
    return result


def bell_coordinates(tensor):
    a,b,c,d = (tensor[q] for q in PATTERNS)
    return sp.Matrix([a+b,a-b,c+d,c-d])


def bell_gradients():
    x = sp.symbols('x:4')
    tensor = bell(x)
    for start in (0,1):
        m = cut(tensor,start)
        adj = m.adjugate().T
        vec = {positions(i,j,start):adj[i,j] for i,j in product(range(4),repeat=2)}
        assert all(sp.expand(a-m.det().diff(v))==0 for a,v in zip(bell_coordinates(vec),x))
    a,b,c,d = x
    assert sp.factor(cut(tensor,0).det()) == a*b*c*d
    assert sp.expand(16*cut(tensor,1).det()+
                     (a-b-c-d)*(a-b+c+d)*(a+b-c+d)*(a+b+c-d)) == 0


if __name__ == '__main__':
    schur_transport()
    bell_gradients()
    print('Both reciprocal transports and the two gradient formulas agree with the actual cuts.')
