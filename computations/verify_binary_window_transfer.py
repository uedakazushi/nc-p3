"""Verify rank-one Hankel graphs and a four-edge chain with nontrivial monodromy."""

from itertools import combinations, combinations_with_replacement, product

import sympy as sp


def hankel(i, j):
    return sp.Matrix(i + 1, j + 1, lambda r, c: values[r + c])


values = sp.symbols("h:6")


def power_vector(d, s, t):
    return sp.Matrix([s ** (d - i) * t**i for i in range(d + 1)])


def check():
    s, t = sp.symbols("s t")
    h = hankel(2, 3)
    points = dict(zip(values, power_vector(5, s, t)))
    assert h.subs(points) == power_vector(2, s, t) * power_vector(3, s, t).T
    minors = [
        sp.expand(h.extract(r, c).det())
        for r, c in product(combinations(range(3), 2), combinations(range(4), 2))
    ]
    standard = [
        values[i] * values[j + 1] - values[i + 1] * values[j]
        for i, j in combinations(range(5), 2)
    ]
    mons = list(combinations_with_replacement(range(6), 2))

    def coefficients(fs):
        return sp.Matrix(
            [
                [
                    sp.Poly(f, *values).coeff_monomial(values[i] * values[j])
                    for i, j in mons
                ]
                for f in fs
            ]
        )

    a, b = coefficients(minors), coefficients(standard)
    assert a.rank() == b.rank() == a.col_join(b).rank() == 10
    shifts = [
        sp.Matrix([[1, 1], [0, 1]]),
        sp.Matrix([[2, 0], [0, 1]]),
        sp.eye(2),
        sp.eye(2),
    ]
    param = sp.Matrix([s, t])
    for i, g in enumerate(shifts):
        degree = 2 + i % 2
        nxt = g * param
        # Each rank-one pair traces a quintic spanning a six-plane.
        forms = list(power_vector(degree, *param) * power_vector(5 - degree, *nxt).T)
        coeff = sp.Matrix(
            [
                [sp.Poly(f, s, t).coeff_monomial(s ** (5 - j) * t**j) for j in range(6)]
                for f in forms
            ]
        )
        assert coeff.rank() == 6
        param = nxt
    assert sp.expand(param[0] * t - param[1] * s) != 0
    print(
        "The catalecticant minors cut out the rational normal quintic; all four edge spans match with nonidentity return on the parameter line."
    )


if __name__ == "__main__":
    check()
