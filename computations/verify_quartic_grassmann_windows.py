"""Dimension and lattice checks for companion quartic maps and elliptic octics."""
from math import comb

import sympy as sp


def main():
    # The numbers are independent of local coordinates of the elliptic curve.
    assert [8*n-comb(n+3, 3) for n in (1, 2, 3)] == [4, 6, 4]
    assert comb(7, 3)-32 == 3
    assert 1+4*(8-4)-3 == 14
    h = sp.symbols('h')
    bundle = 2-4*h**2
    exp2 = sp.exp(2*h).series(h, 0, 4).removeO()
    exp_minus = sp.exp(-2*h).series(h, 0, 4).removeO()
    assert sp.rem(sp.expand(bundle*(exp2+exp_minus)-4), h**4, h) == 0
    # In the quartic surface plane net C and C' are residual (4,4) curves.
    gram = sp.Matrix([[4, 8], [8, 0]])
    polarization, curve = sp.eye(2)[:, 0], sp.eye(2)[:, 1]
    residual = 4*polarization-curve
    assert (polarization.T*gram*residual)[0] == 8
    assert (residual.T*gram*residual)[0] == 0
    assert (curve.T*gram*residual)[0] == 32
    print('PASS: elliptic octic Rao dimensions, companion characters and residual lattice')


if __name__ == '__main__':
    main()
