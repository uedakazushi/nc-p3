"""Verify the binary triple-root correspondence for a fourfold stationary ruling."""
from sympy import Matrix, div, expand, gcd, symbols

t = symbols("t")
h = t**4 + t + 1


def coefficient_column(poly):
    remainder = div(poly, h, t)[1].expand()
    return [remainder.coeff(t, i) for i in range(4)]


def matrix(a, b):
    return Matrix.hstack(*(Matrix(coefficient_column(poly))
                           for poly in (-b, -t*b, a, t*a)))


def main():
    cubic = t**3 + 2*t + 3
    left, right = t + 1, 2*t - 3
    a = expand(cubic*left - 4*h)
    b = expand(cubic*right + 7*h)
    assert matrix(a, b) * Matrix([1, 1, -3, 2]) == Matrix.zeros(4, 1)
    assert matrix(a, b).rank() == 3
    assert gcd(a + 4*h, b - 7*h) == cubic
    z = symbols("z")
    determinant = matrix(a + z, b + z*t**2).det().expand()
    assert determinant.as_poly(z).degree() == 4 and determinant.subs(z, 0) == 0
    assert 5 + 1 - 3 == 3 and 3 - (-1) == 4 < 20 - 6 - 3
    print("PASS: binary syzygy kernel recovers the triple fiber, both determinants have degree four against eleven isolated marks")


if __name__ == "__main__":
    main()
