"""Check the common central degree and the path bounds for the domain branch."""
from math import lcm

import sympy as sp


def main():
    degree = 4*lcm(1,2,3,4)
    socle = 4*(degree-1)
    assert (degree, socle+1, 3*max(socle, degree), socle+4) == (48,189,564,192)
    t = sp.symbols("t")
    coefficients = (sp.Poly(sum(t**j for j in range(degree)), t)**4).all_coeffs()
    assert len(coefficients) == 189 and coefficients == coefficients[::-1]
    assert sum(coefficients) == degree**4
    print("Four characteristic coefficients, degree-48 parameters, cutoff degrees 189, 564, and 192.")


if __name__ == "__main__":
    main()
