"""Section Euler characteristics, character parity, and the A1 multiplicity budget."""
from itertools import product
from sympy import Rational, expand, symbols


def main():
    m, a = symbols("m a")
    x = 4*m + a
    hilbert = (x+3)*(x+2)*(x+1)*Rational(1, 6)
    cut = expand(hilbert - hilbert.subs(m, m-2))
    assert cut == 64*m**2 + 32*a*m - 64*m + 4*a**2 - 16*a + 36
    assert cut.subs({a:0, m:0}) == 1+35
    assert [signs for signs in product((-1, 1), repeat=2)
            if 1 not in signs and signs[0]*signs[1] == -1] == []
    types = [(d, r) for d in (2, 3, 4) for r in (1, 2) if d*r <= 4]
    assert [(x, y) for x, y in product(types, repeat=2)
            if x[0]*x[1]+y[0]*y[1] == 4] == [((2, 1), (2, 1))]
    print("Surface Euler polynomial, two-dimensional determinant character, and balanced A1 multiplicities checked.")


if __name__ == "__main__":
    main()
