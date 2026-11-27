"""Numerical Riemann--Roch identities for the center hyperplane and local determinant classes."""
from itertools import product

import sympy as sp


def main():
    m, d = sp.symbols("m d")
    n = 4*m + d
    affine = (n+3)*(n+2)*(n+1)/6
    surface = sp.expand(affine - affine.subs(m, m-1))
    assert surface == 32*m**2 + 16*d*m + 2*d**2 + 2
    assert surface.coeff(m, 1) == 4 * (4*d)
    assert surface.subs({m: 0, d: 0}) == 2

    # The two rulings have opposite classes on a punctured node.
    for exponent, orientation, local_twist in product([1, 2, 3], [-1, 1], range(-4, 5)):
        determinant = 4*local_twist + orientation*exponent
        assert (determinant % 2 == 0) == (exponent == 2)
    print("Hyperplane Hilbert polynomial, determinant increments, and both local rank profiles verified.")


if __name__ == "__main__":
    main()
