"""Conductor values and first jets of the normalized residual identifications."""

import sympy as sp


def main():
    r,s,z,c,t,lam = sp.symbols('r s z c t lam', nonzero=True)
    node = (r*z-c*s)/(r*z-s)
    assert node.subs(s,0) == 1 and sp.cancel(node.subs(r,0)) == c
    branch = lam-s/(r*z-s)
    assert branch.subs(s,0) == lam and branch.subs(r,0) == lam+1
    cusp = (r-(1+2*t*r)*z)/(r-z)
    assert cusp.subs(r,0) == 1 and sp.simplify(sp.diff(cusp,r).subs(r,0)) == 2*t
    assert sp.simplify(sp.diff(cusp**2,r).subs(r,0)) == 4*t
    print('Verified the two nodal conductor multipliers, the additive shift, and the cuspidal first jets.')


if __name__ == '__main__':
    main()
