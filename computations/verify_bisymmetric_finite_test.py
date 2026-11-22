"""The reflected binary quotient equations and the two uniform path bounds."""

from itertools import product
from math import comb

import sympy as sp

from verify_conductor_pencil import coefficients, t


def opposite_bases(relations, rows):
    """Linear constraints on the four column vectors after choosing three independent rows."""
    assert rows.shape == (3, 4) and rows.rank() == 3
    equations = []
    for relation in relations:
        restriction = rows * relation
        for i, j in product(range(2), range(3)):
            row = sp.zeros(1, 16)
            row[:, 4 * j : 4 * j + 4] = restriction[i + 1, :]
            row[:, 4 * (j + 1) : 4 * (j + 1) + 4] = -restriction[i, :]
            equations.append(row)
    return sp.Matrix.vstack(*equations)


def check():
    pencil = [1 + 2 * t + 4 * t**2 + 7 * t**6 + 3 * t**7, t + 4 * t**4 + t**7]
    quotient = sp.Matrix.hstack(*coefficients(pencil, 7).T.nullspace()).T
    cubics = [1, t, t * t, t**3]
    rows = sp.eye(4)[1:, :]
    basis = sp.Matrix([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1], [1, 0, 0, 2]])
    assert basis.det() != 0
    for quartics in ([1 + t**4, t, t * t, t**3], [1, t * t, t**3, t**4]):
        forms = quotient * coefficients(
            [f * g for f, g in product(quartics, cubics)], 7
        )
        relations = [forms.row(i).reshape(4, 4) * basis for i in range(6)]
        equations = opposite_bases(relations, rows)
        expected = basis.inv().T.reshape(16, 1)
        assert equations.rank() == 15 and equations * expected == sp.zeros(36, 1)
    for n, standard, rectangle, binary in [
        (17, 1140, 512, 1152),
        (21, 2024, 2048, 4608),
    ]:
        assert comb(n + 3, 3) == standard
        assert 2 ** ((n + 1) // 2) == rectangle
        assert 9 * 2 ** ((n - 3) // 2) == binary > standard
    assert 64 - 2 * 24 + 4 == comb(6, 3) == 20
    print(
        "The two conductor models each recover a unique invertible column basis up to scale. Degree 21 detects both rectangular and binary quotients."
    )


if __name__ == "__main__":
    check()
