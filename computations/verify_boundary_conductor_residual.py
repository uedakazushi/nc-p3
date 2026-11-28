"""Polynomial identities for the residual of the conductor in singular boundary windows."""

import sympy as sp

from verify_singular_boundary_bezout import divide


def check(kind):
    z, w, u = sp.symbols('z w u')
    a0, a1, a2, a3 = sp.symbols('a0 a1 a2 a3')
    if kind == 'node':
        a3 = u*a0
    else:
        a1 = u*a0
    az = a0+a1*z+a2*z**2+a3*z**3
    aw = az.subs(z, w)
    if kind == 'node':
        uz = sp.cancel(u*(az-a0)/z-z*(az-a3*z**3)/u)
        particular = (z+w)*az*aw+divide(a3*(z*az-w*aw)+a0*z*w*(z**4*aw-w**4*az), z, w)
        square = az*aw*(u+z)*(u+w)/u
        conductor = z*w
    else:
        uz = sp.cancel((a1+(-u+u**2*z)*az)/z**2)
        particular = az*aw-a1*(z*az+w*aw)-a1*z*w*divide(az-aw, z, w)
        square = az*aw*(1-u*z)*(1-u*w)
        conductor = (z*w)**2
    bezout = divide(az*uz.subs(z, w)-aw*uz, z, w)
    assert sp.cancel(particular+conductor*bezout-square) == 0

    c0, c1, c2 = sp.symbols('c0 c1 c2')
    qz = c0+c1*z+c2*z**2
    g = z if kind == 'node' else z**2
    small = divide(g*qz.subs(z, w)-g.subs(z, w)*qz, z, w)
    expected = c0-c2*z*w if kind == 'node' else c0*(z+w)+c1*z*w
    assert sp.expand(small-expected) == 0


if __name__ == '__main__':
    check('node')
    check('cusp')
    print('Verified both universal Schur identities and both conductor Bezout subspaces.')
