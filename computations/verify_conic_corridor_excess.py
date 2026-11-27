"""Exact coefficient checks for the conic corridor excess criterion."""
from itertools import product

from sympy import Matrix, symbols


def main():
    s, t = symbols("s t")
    # A representative L2 plus a regular block. Its last row is a pure
    # relation, and its last letter is the plane equation of the kernel conic.
    pencil = Matrix([[s, t, 0, 0], [0, s, t, 0], [0, 0, 0, s + t]])
    kernel = Matrix([t**2, -s*t, s**2, 0])
    assert pencil * kernel == Matrix.zeros(3, 1)
    assert pencil[2, :] == (s + t) * Matrix([[0, 0, 0, 1]])
    assert Matrix([[value.coeff(s, a).coeff(t, 2-a) for value in kernel]
                   for a in range(3)]).rank() == 3
    restricted = Matrix.hstack(*(Matrix.vstack(-b*pencil.subs({s:a, t:b}),
                                                a*pencil.subs({s:a, t:b}))
                                 for a, b in [(1, 0), (0, 1)]))
    assert restricted.rank() == 6
    first = Matrix(3, 4, symbols("a:12"))
    second = Matrix(3, 4, symbols("b:12"))
    xi = Matrix(1, 3, symbols("x:3"))
    row = xi * (s*first + t*second)
    tails = Matrix.vstack(*[factor.applyfunc(lambda v: v.expand().coeff(s, a).coeff(t, 2-a))
                            for factor, a in product((s*row, t*row), range(3))])
    assert tails == Matrix.vstack(Matrix.zeros(1, 4), xi*second, xi*first,
                                  xi*second, xi*first, Matrix.zeros(1, 4))
    print("PASS: conic plane, six-dimensional restriction and two tail coefficients")


if __name__ == "__main__":
    main()
