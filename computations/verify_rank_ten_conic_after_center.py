"""Check the left relations and the projected point pairs of both pencils."""
import sympy as sp


def restriction_rank(matrix, s, t):
    columns = []
    for col in [matrix.col(j) for j in range(matrix.cols)]:
        stacked = sp.Matrix.vstack(s*col, t*col)
        for monomial in (s*s, s*t, t*t):
            columns.append(stacked.applyfunc(lambda x: sp.Poly(x, s, t).coeff_monomial(monomial)))
    return sp.Matrix.hstack(*columns).rank()


def main():
    s, t, a, b, c = sp.symbols("s t a b c")
    binary = sp.Matrix([[s, 0], [t, s], [0, t]])
    ternary = sp.Matrix([[s, t, 0], [0, 0, s], [0, 0, t]])
    assert sp.Matrix([[t*t, -s*t, s*s]]) * binary == sp.zeros(1, 2)
    assert restriction_rank(binary, s, t) == 6
    assert restriction_rank(ternary, s, t) == 5
    assert ternary.rank() == 2
    assert ternary * sp.Matrix([-t, s, 0]) == sp.zeros(3, 1)
    expression = (sp.Matrix([[c, -a, -b]]) * ternary * sp.Matrix([a, b, c]))[0]
    assert sp.expand(expression) == 0
    print("PASS: binary exclusion, singular ternary pencil and both planar projections")


if __name__ == "__main__":
    main()
