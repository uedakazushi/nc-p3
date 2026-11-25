"""Check intersection, Hilbert-polynomial and determinant identities exactly."""
from math import comb

from sympy import Poly, Rational, cyclotomic_poly, expand, symbols


def projective_three_euler(degree):
    return (degree + 1) * (degree + 2) * (degree + 3) * Rational(1, 6)


def main():
    a, b = symbols("a b")
    hilbert = sum(
        (-1) ** j * comb(5, j)
        * projective_three_euler(a - j) * projective_three_euler(b - j)
        for j in range(6)
    )
    assert expand(hilbert) == 10 * a + 10 * b - 10
    determinant = [sum((-1) ** j * comb(4, j) for j in range(i + 1))
                   for i in range(4)]
    assert determinant == [1, -3, 3, -1]
    assert [Poly(cyclotomic_poly(m, a), a).eval(1) for m in (2, 3, 4, 6)] == [2, 3, 2, 1]
    assert comb(6, 3) - (4 * 4 - 2 * 1 + 2) == 4
    print("PASS: (10,10), genus 11, third difference, torsion orders and excess length")


if __name__ == "__main__":
    main()
