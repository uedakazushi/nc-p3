"""The three exterior identities and their nonzero Pfaffian multiples."""

from itertools import combinations

import sympy as sp


def exterior(u, v):
    return u*v.T-v*u.T


def check(kind):
    z, w = sp.symbols('z w')
    t = sp.Rational(3, 2)
    a = 1+2*z-3*z**2+4*z**3
    b = -2+3*z+5*z**2-2*z**3
    if kind == 'node':
        a += t*z**4
        b += -2/t*z**4
        multiplied = (z-t*w)*a*b.subs(z, w)+(t*z-w)*b*a.subs(z, w)
    elif kind == 'cusp':
        a += (t-2)*z+3*z**4
        b += (2*t-3)*z-2*z**4
        multiplied = (z-w)*(a*b.subs(z, w)+b*a.subs(z, w))/t+z*w*(a*b.subs(z, w)-b*a.subs(z, w))
    else:
        a -= z**4
        b += 3*z**4
        multiplied = (z-w)*a*a.subs(z, w)+2*(z+w)*(b*a.subs(z, w)-a*b.subs(z, w))
    columns = [sp.Matrix([sp.expand(f).coeff(z, i) for i in range(6)]) for f in [z*a, a, z*b, b]]
    v, a0, q, b0 = columns
    if kind == 'node':
        skew = exterior(v, b0)+t*exterior(q, a0)
    elif kind == 'cusp':
        skew = (exterior(v, b0)+exterior(q, a0))/t+exterior(v, q)
    else:
        skew = exterior(v, a0)+2*exterior(q, a0)-2*exterior(v, b0)
    matrix = sp.Poly(multiplied, z, w)
    assert skew == sp.Matrix(6, 6, lambda i, j: matrix.coeff_monomial(z**i*w**j))
    pencil = sp.Matrix.hstack(*columns)
    factor = None
    for indices in combinations(range(6), 4):
        i, j, k, l = indices
        pf = skew[i,j]*skew[k,l]-skew[i,k]*skew[j,l]+skew[i,l]*skew[j,k]
        minor = pencil.extract(indices, range(4)).det()
        if minor and factor is None:
            factor = pf/minor
        assert pf == (factor or 0)*minor
    assert factor
    return factor


if __name__ == '__main__':
    print({kind: check(kind) for kind in ('node', 'cusp', 'branch')})
