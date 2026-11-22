"""Certify the four formal rank-fifteen classes outside the ordinary multiplier chart."""

import sympy as sp
from flint import fq_default_ctx

from verify_klein_center_lift import echelon, field_paths
from verify_rank_fifteen_parameter_orbits import amplitudes, pattern, product_sign
from verify_rank_fifteen_pencil_linear_residues import graded_profile
from verify_rank_fifteen_regular_orbits import tensor_windows
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile


def dual_relations(windows, unit):
    dual = []
    for side in windows:
        relations = []
        for h in range(4):
            pivots = echelon(
                [
                    {i: c for (i, j), c in row.items()}
                    for color, row in side
                    if color == h
                ]
            )
            for j in range(4):
                if j in pivots:
                    continue
                vector = {(j, j ^ h): unit}
                for i, row in pivots.items():
                    if j in row:
                        vector[i, i ^ h] = -row[j]
                relations.append((h, vector))
        assert len(relations) == 10
        dual.append(relations)
    return dual


def verify_residuals():
    field = fq_default_ctx(17, 4, "v")
    one = field.one()
    unit = Scalar(one)
    root = field(2).sqrt()
    for K in (field(0), field(3), one + root, one - root):
        U = K * K
        H = (U - 1) / 2 if K in (field(0), field(3)) else -K
        V = H * H
        x = (U + 3) / 4
        assert (U - 1) ** 2 == 4 * V and V * (U + 3) / (U - 1) ** 2 == x
        r = K + H + ((K + H) ** 2 - 1).sqrt()
        w = K - H + ((K - H) ** 2 - 1).sqrt()
        s, t, n = r.sqrt(), w.sqrt(), ((U - 1) / 2).sqrt()
        gamma = amplitudes(s, t / s, n, one)
        assert all(not c.is_zero() for c in gamma[1:])
        windows = tensor_windows(gamma, field)
        basis, _ = field_paths(dual_relations(windows, unit), unit)
        for p in range(4):
            assert [sum(len(basis[p, j, h]) for h in range(4)) for j in range(6)] == [
                1,
                4,
                6,
                4,
                1,
                0,
            ]
        try:
            assert graded_profile(windows, unit)
        except AssertionError:
            assert verify_central_profile(windows, unit)


def verify_identities():
    K, H, x, y, v = sp.symbols("K H x y v")
    cubic = v * (x - 4 * v) ** 2 - y
    assert sp.factor(sp.discriminant(cubic, v) - 256 * y * (x**3 - 27 * y)) == 0
    assert (
        sp.factor(
            ((K + H) ** 2 - 1) * ((K - H) ** 2 - 1)
            - ((K * K + H * H - 1) ** 2 - 4 * K * K * H * H)
        )
        == 0
    )
    boundary = (((K + H) ** 2 - 1) * ((K - H) ** 2 - 1)).subs(H, (K * K - 1) / 2)
    assert sp.factor(boundary) == (K - 3) * (K - 1) ** 3 * (K + 1) ** 3 * (K + 3) / 16
    a, b, L, M, A, B, C, D = sp.symbols("a b L M A B C D")
    src = pattern(a, b, L, M, A, B, C, D)
    # Reversing a word transposes both factors, exchanges them and reindexes their phases.
    reversed_row = [
        (-1) ** product_sign(j, j) * src[j % 4 + 4 * (0, 3, 2, 1)[j // 4]]
        for j in range(16)
    ]
    assert tuple(reversed_row) == pattern(a, b, L, M, A, B, -D, -C)


if __name__ == "__main__":
    verify_identities()
    verify_residuals()
    print(
        "Four excluded quotient points: dual 14641 and four-start central certificates at all four"
    )
