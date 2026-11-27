"""Check the free blocks and all shuffle signs in the Klein-fixed bisymmetric locus."""
from itertools import product

import sympy as s

K = range(4)


def orbit(i,j,k,n):
    return {(i,j,k,n),(i,n,k,j),(k,j,i,n),(k,n,i,j)}


def main():
    a,b = s.symbols('a:4'),s.symbols('b:4')
    values = {}
    for i,k in product(K,repeat=2):
        for q in orbit(i,i^2,k,k^3):
            assert q not in values
            values[q] = a[i]*b[k]
    # Some thin entries have already specified the off-diagonal thick blocks.
    free = []
    for i,k in product(K,repeat=2):
        q = (i,i,k,k^1)
        if q not in values:
            v = s.Symbol(f'x{i}{k}')
            free.append(v)
            for u in orbit(*q):
                assert u not in values
                values[u] = v
    assert len(free)==8 and len(values)==64
    blocks = [s.Matrix(4,4,lambda i,k:values[i,i^h,k,k^h^1]) for h in K]
    assert blocks[1]==blocks[0].T and blocks[3]==blocks[2].T
    U,V = blocks[0][:2,2:],blocks[0][2:,:2]
    assert U==s.Matrix([[a[2]*b[0],a[0]*b[3]],[a[1]*b[2],a[3]*b[1]]])
    assert V==s.Matrix([[a[0]*b[2],a[2]*b[1]],[a[3]*b[0],a[1]*b[3]]])
    assert U.det()==-V.det() and U.det()!=0
    def sign(i,j):
        return -1 if (i>>1)&(j&1) else 1
    def twist(q):
        return s.prod(sign(q[i],q[j]) for i in K for j in range(i+1,4))
    for q in product(K,repeat=4):
        i,j,k,n = q
        h = i^j^k^n
        comm = sign(i^j,h)*sign(h,i^j)
        assert twist(q)*twist((k,n,i,j))==comm
        assert sign(i,j)*sign(j,i)==sign(i,i)*sign(j,j)*sign(i^j,i^j)
    print('PASS: twenty-four independent coordinates, both determinant blocks and both cocycle identities')


if __name__=='__main__':
    main()
