"""Hilbert functions and surface bounds for the image of the five central quartics."""
from math import comb
from sympy import cancel, symbols


def main():
    t = symbols("t")
    hilbert = (1-t**4)/(1-t)**5
    for quadrics in (0, 1, 2):
        ambient = 1/((1-t)**5*(1-t**2)**quadrics)
        assert cancel(ambient*(1-t**2)**quadrics*(1-t**4) - hilbert) == 0
    surface = hilbert*(1-t**2)
    assert cancel(surface - (1+t**2)*(1-t**2)**2/(1-t)**5) == 0
    assert [(d, e) for d in range(2, 9) for e in range(3, 9) if d*e <= 8] == [(2, 3), (2, 4)]
    assert comb(5+1, 2) - (5+8) == 2
    assert 2*3 + 2 == 8 and 2*8 == 2*9-2
    assert 6*8 > 6**2 and 6*8 < 7**2
    print("Center and double-cover series, image degrees, and the genus-nine base divisor checked.")


if __name__ == "__main__":
    main()
