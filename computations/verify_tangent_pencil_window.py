"""Two finite symmetric windows on opposite sides of the periodic PBW condition."""

from itertools import combinations

import sympy as sp


def rank_minors(matrix):
    indices = list(combinations(range(4), 3))
    return [
        matrix.extract(i, j).det() for m, i in enumerate(indices) for j in indices[m:]
    ]


def verify():
    s, t, u, v = sp.symbols("s t u v")
    q0 = sp.diag(1, 1, 0, 0)
    q1 = sp.zeros(4)
    q1[:2, 2:] = q1[2:, :2] = sp.eye(2)
    q2 = sp.Matrix([[2, 3, 1, -1], [3, -1, 0, 2], [1, 0, 1, 0], [-1, 2, 0, -1]])
    q3 = sp.Matrix([[0, 2, -1, 3], [2, 1, 2, 1], [-1, 2, 0, 1], [3, 1, 1, 0]])
    matrices = [q0, q1, q2, q3]
    letters = sp.symbols("a:d")
    vector = sp.Matrix(letters)
    equations = [(vector.T * q * vector)[0] for q in matrices]
    assert sp.groebner(
        equations, *letters, domain=sp.QQ, order="grevlex"
    ).is_zero_dimensional
    matrix = s * q0 + t * q1 + u * q2 + v * q3
    minors = rank_minors(matrix)
    infinity = sp.groebner(
        [f.subs(s, 0) for f in minors], t, u, v, domain=sp.QQ, order="grevlex"
    )
    assert infinity.is_zero_dimensional
    # Make t+3u+7v the last (separating) coordinate of the affine algebra.
    z = sp.Symbol("z")
    affine = [f.subs({s: 1, t: z - 3 * u - 7 * v}) for f in minors]
    lex = sp.groebner(affine, u, v, z, domain=sp.QQ, order="grevlex").fglm("lex")
    assert [p.LM(order="lex").exponents for p in lex.polys] == [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 10),
    ]
    polynomial = sp.Poly(lex[-1], z)
    assert polynomial.gcd(polynomial.diff()).monic().as_expr() == z
    assert polynomial.nth(0) == polynomial.nth(1) == 0 and polynomial.nth(2) != 0
    # A rank-one point would have zero gradients for every 3x3 minor.
    # All points except the origin above are simple. At the origin q0 has rank two.
    print(
        "PASS basepoint-free web; rank-two section has eight simple points and one double point"
    )
    print("PASS ordered edge has sixteen simple points and two double points")


def triangular_window():
    x, y, z, w = letters = sp.symbols("x y z w")
    forms = [
        x * x + 2 * x * y + 6 * x * z + 2 * y * z + 4 * x * w - 2 * y * w + 2 * z * w,
        y * y + 2 * x * z + 4 * y * z + 6 * x * w + 2 * y * w - 2 * z * w,
        z * z + 2 * x * w + 6 * y * w + 4 * z * w,
        w * w,
    ]
    matrices = [sp.hessian(f, letters) / 2 for f in forms]
    a, b, c, d = sp.symbols("a b c d")
    matrix = sum(
        (t * m for t, m in zip((a, b, c, d), matrices, strict=True)), sp.zeros(4)
    )
    minors = rank_minors(matrix)
    affine = sp.groebner(
        [f.subs(a, 1) for f in minors], b, c, d, domain=sp.QQ, order="grevlex"
    ).fglm("lex")
    assert [p.LM(order="lex").exponents for p in affine.polys] == [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 5),
    ]
    p = sp.Poly(affine[-1], d)
    assert p.gcd(p.diff()).degree() == 0
    for chart, variables in (({a: 0, b: 1}, (c, d)), ({a: 0, b: 0, c: 1}, (d,))):
        assert list(
            sp.groebner([f.subs(chart) for f in minors], *variables, domain=sp.QQ)
        ) == [1]
    assert matrices[-1].rank() == 1
    print(
        "PASS triangular PBW web; five reduced rank-two quadrics and one rank-one support"
    )
    print("PASS ordered edge has ten reduced points and one point of length ten")


if __name__ == "__main__":
    verify()
    triangular_window()
