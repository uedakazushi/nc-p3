"""Three local-GL contractions recover the full-rank transverse parameters."""
from itertools import product

import numpy as np
import sympy as s

from verify_klein_three_lines import K, twist


def universal(a,b,c,d):
    A,C = s.Matrix([a,b,c,d]),s.Matrix([c,d,a,b])
    first = s.Matrix.hstack(A,A,C,C)
    mixed = s.Matrix([c,d,c,d])*s.ones(1,4)
    blocks = (first,first.T,mixed,mixed.T)
    w = {}
    for h,i,k in product(K,repeat=3):
        key = (i,i^h,k,k^h^1)
        w[key]=twist(key)*blocks[h][i,k]
    return w


def opposite(w,b):
    return s.Matrix([[w[i,j,i^b,j^b^1] for j in K] for i in K])


def adjacent(w,b):
    return s.Matrix([[w[i,i^b,k,k^b^1] for k in K] for i in K])


def cubic_contraction(w,inverse,p,q,target):
    variables = (1,p,q)
    data = []
    for tensor in (w,{key:s.factor(4*(p+q)*z) for key,z in inverse.items()}):
        coefficients = [np.zeros((4,)*4,dtype=np.int64) for _ in variables]
        for key,val in tensor.items():
            poly = s.Poly(val,p,q)
            assert poly.total_degree()<=1
            for array,mon in zip(coefficients,variables):
                array[key]=int(poly.coeff_monomial(mon))
        data.append(coefficients)
    # After clearing the cubic denominator the bidegree is at most (6,6).
    # Exact integer interpolation by evaluation on a unisolvent 7 x 7 grid.
    numerator = s.cancel(target*(4*(p+q))**3)
    assert s.Poly(numerator,p,q).total_degree()<=6
    for i,j in product(range(7),repeat=2):
        W,V = [a+i*b+j*c for a,b,c in data]
        upper = 4**12*max(1,int(np.max(abs(W))))**3*max(1,int(np.max(abs(V))))**3
        assert upper<2**63
        val = np.einsum('abcd,efgh,ijkl,abgh,ejkd,ifcl->',W,W,W,V,V,V,optimize='greedy')
        assert int(val)==2*numerator.subs({p:i,q:j})


def main():
    p,q,z = s.symbols('p q z')
    t,e = p+q,p-q
    w = universal(p,q,1,1)
    inverse = {}
    for b in K:
        m = opposite(w,b).inv().T.applyfunc(s.factor)
        for i,j in product(K,repeat=2):
            inverse[i,j,i^b,j^b^1]=m[i,j]
    blocks = [(adjacent(w,b)*adjacent(inverse,b).T).applyfunc(s.factor) for b in K]
    A = (t+4/t)**2
    for b,B in enumerate(blocks):
        expected = z**2*((z-2)**2-A/4) if b<2 else z**3*(z-4)
        assert s.factor(B.charpoly(z).as_expr()-expected)==0
    values = [s.factor(sum((B**2)[i,i] for B in blocks)) for i in K]
    assert values[::2]==[values[0]]*2 and values[1::2]==[values[1]]*2
    assert s.factor(values[0]-values[1]-e*(t-4/t)/4)==0
    C = -e**2*(1+4/t**2)+s.Rational(7,2)*(t**2+16/t**2)+100
    cubic_contraction(w,inverse,p,q,C)
    neg = universal(-p,-q,1,1)
    for key,val in w.items():
        assert -s.I**sum(i//2 for i in key)*val==neg[key]
    print('PASS: characteristic polynomial, reduced partial trace, cubic contraction and simultaneous sign')


if __name__=='__main__':
    main()
