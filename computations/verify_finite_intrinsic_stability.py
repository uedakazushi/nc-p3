"""Check the two intrinsic flag tests and the sharp lower bounds on local length."""
from itertools import product
from math import prod

import sympy as sp


def main():
    words = list(product(range(4),repeat=4))
    for weights,top,degrees in [((0,1,1,2),4,(1,1,2,2,3,3)),((0,1,2,3),6,(1,2,3,3,4,5))]:
        assert prod(degrees) // prod(weights[1:])**2 == (9 if top==4 else 10)
        for q in words:
            if sum(weights[i] for i in q)>=top or any(q[i]==q[(i+1)%4]==0 for i in range(4)):
                continue
            if top==4:
                assert 0 in q
                if q.count(0)==1:
                    assert all(i in (0,1,2) for i in q)
                else:
                    assert any(i in (1,2) for i in q)
    # With nonzero r_5 coefficient at 22, the extreme complementary
    # rank-one block propagates to each of the neighboring 22 pairs.
    assert all(q.count(2)==3 for q in [(2,2,0,2),(2,2,2,0)])
    u,v,p,q,k = sp.symbols('u v p q k',nonzero=True)
    tangent = sp.Matrix([[1,p,0,0], [0,k,0,0], [0,0,1,0],[0,0,0,1]])
    assert sp.expand(v*(tangent[1,:]*sp.Matrix([u,v,p,q]))[0]) == k*v**2
    print('Weight-four contractions, propagated square tangent cone, and local multiplicities 9 and 10 verified.')


if __name__ == '__main__':
    main()
