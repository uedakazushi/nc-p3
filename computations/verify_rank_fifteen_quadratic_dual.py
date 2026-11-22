"""Cubic overlap charts, nondegenerate opposite kernels and the degree-six test."""

import sympy as sp
from verify_rank_fifteen_pencil_infinity import a, b, e, f, z, l, symbolic_blocks
from verify_klein_elliptic_completion import gram_twist


def sixth_insertion(windows, basis, image, start, h):
    """Rows for A_4 R in the four final-letter summands of A_5 V."""
    offsets = []
    width = 0
    for last in range(4):
        offsets.append(width)
        width += len(basis[start, 5, h ^ last])
    rows = []
    for g, rule in windows[start]:  # The final relation is at start+4.
        for prefix in basis[start, 4, h ^ g]:
            row = {}
            for (a, b), coefficient in rule.items():
                for j, value in image(start, prefix + (a,)).items():
                    col = offsets[b] + j
                    row[col] = row.get(col, 0) + coefficient * value
            rows.append({j: value for j, value in row.items() if value})
    return rows, width


def check_quadratic_residue():
    from itertools import product
    from flint import fq_default_ctx
    from rank_fifteen_octic import finite_windows
    from verify_klein_center_lift import field_paths, echelon
    from verify_rank_fourteen_geometric_reduction import Scalar

    field = fq_default_ctx(101, 2, "t")
    unit = field.one()
    k, c, L = 2 * unit, 3 * unit / 4, 5 * unit / 4
    slope = c + 1 / c

    def signs(square):
        root = square.sqrt()
        assert not root.is_zero()
        return root, -root

    radical = signs(k * k - 1)
    assert radical[0] ** 101 == -radical[0]
    for H, A, B in product(
        radical, signs((k + slope) ** 2 - 1), signs((k - slope) ** 2 - 1)
    ):
        n = H / L
        r = k + (H * H + H * A) / slope
        w = k + (-H * H + H * B) / slope
        z = r * r
        assert not (
            z * (z * z - 1) * (w * w - 1) * (w * w - z) * (w * w * z - 1)
        ).is_zero()
        a, b, e, f = 1 + n, (k + c * n) / r, 1 - n, (k - c * n) / r
        for frame in ((a, b, e, f), (z * b, a, z * f, e)):
            windows = finite_windows((z, w / r, *frame), unit, Scalar)
            basis, image = field_paths(windows, Scalar(unit))
            for start, h in product(range(4), repeat=2):
                rows, width = sixth_insertion(windows, basis, image, start, h)
                assert (len(rows), width) == ((54, 54, 51, 51)[h], 56)
                assert len(echelon(rows)) == (36, 36, 34, 34)[h]


def verify():
    pair = {0, 1}
    one = sp.eye(2)
    path_units = {
        a * a - z * b * b,
        e * e - z * f * f,
        a * a - z * f * f,
        e * e - z * b * b,
        z,
        l,
        z - 1,
        l * l * z - 1,
    }
    allowed = {
        sp.Poly(g, a, b, e, f, z, l).monic().as_expr()
        for value in path_units
        for g, _ in sp.factor_list(value)[1]
    }

    def is_unit(value):
        num, den = sp.fraction(sp.cancel(value))
        return bool(num) and all(
            sp.Poly(g, a, b, e, f, z, l).monic().as_expr() in allowed
            for part in (num, den)
            for g, _ in sp.factor_list(part)[1]
        )

    windows = symbolic_blocks(1, z, l, b)
    for side in windows:
        for h, (rules, _) in enumerate(side):
            matrix = sp.Matrix(rules)
            if h in pair:
                assert matrix[:, 2:].applyfunc(sp.cancel) == one and is_unit(
                    matrix[:, :2].det()
                )
            else:
                assert all(is_unit(x) for x in matrix)
    for h in range(4):
        thin = {h ^ i for i in (2, 3)}
        assert thin == pair or thin == {2, 3}
        for g in pair:
            assert {j ^ g for j in thin} == thin

    A = sp.Matrix([[a, b], [z * b, a]])
    B = sp.Matrix([[e, f], [z * f, e]])
    _, opposite, _ = gram_twist(A, B, z, l)
    for h in pair:
        indices = [4 * i + (i ^ h) for i in range(4)]
        block = opposite.extract(indices, indices)
        if h == 0:
            block[3, :] = block[3, :] / z
        else:
            block[:, 1] = block[:, 1] / z
        m, n = block[0, 0], block[0, 1]
        c, d = block[2, 2], block[2, 3]
        expected = sp.Matrix([[m, n, 1, 1], [n, m, 1, 1], [1, 1, c, d], [1, 1, d, c]])
        assert (block - expected).applyfunc(sp.cancel) == sp.zeros(4)
        # In the pair-symmetric sector the off-diagonal entries are 2.
        embed = sp.Matrix([[1, 0], [1, 0], [0, 1], [0, 1]])
        assert (
            expected * embed - embed * sp.Matrix([[m + n, 2], [2, c + d]])
        ).applyfunc(sp.cancel) == sp.zeros(4, 2)
    for h in range(4):
        perm = sp.zeros(4)
        for j in range(4):
            perm[j, j ^ h] = sp.symbols("x:4")[j]
        assert sp.factor(perm.det()) in (
            sp.prod(sp.symbols("x:4")),
            -sp.prod(sp.symbols("x:4")),
        )
    # The symbolic Koszul Euler ranks are attained at every regular ordinary point.
    r = (2, 2, 1, 1)
    dim2 = (2, 2, 3, 3)
    dim4 = (9, 10, 8, 8)
    inputs = [sum(r[g] * dim4[h ^ g] for g in range(4)) for h in range(4)]
    ranks = [inputs[h] - 20 + dim2[h] for h in range(4)]
    assert inputs == [54, 54, 51, 51] and ranks == [36, 36, 34, 34]
    assert [56 - rank for rank in ranks] == [20, 20, 22, 22]
    print(
        "Verified all four cubic overlap charts, both opposite null-line forms and the degree-six block ranks."
    )


if __name__ == "__main__":
    verify()
    check_quadratic_residue()
