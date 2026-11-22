"""Transport endpoint and octic certificates, including nonordinary target parameters."""

from itertools import product

import sympy as sp
from flint import fq_default_ctx

from verify_klein_center_lift import echelon
from verify_rank_fifteen_parameter_orbits import (
    amplitudes,
    dot,
    recover,
    template_group,
)
from verify_rank_fifteen_pencil_linear_residues import graded_profile
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile


def square_roots(value):
    root = value.sqrt()
    return (root,) if root.is_zero() else (root, -root)


def lifted_coordinates(K, h):
    """r,1/r have sum 2(K+h); w,1/w have sum 2(K-h)."""
    for x, y in product(square_roots((K + h) ** 2 - 1), square_roots((K - h) ** 2 - 1)):
        yield from product(square_roots(K + h + x), square_roots(K - h + y))


def exceptional_parameters(field, octic):
    one = field.one()
    if octic:
        levels = [field(0), field(2), field(-2)]
        levels.extend(
            (a + b) / 2 for a, b in product((one, -one), square_roots(field(5)))
        )
        for K in levels:
            for eta in square_roots(field(-3)):
                for n in square_roots((K * K - 1) * (1 + eta) / 2):
                    for h in (K * K - 1, 1 - K * K):
                        for s, t in lifted_coordinates(K, h):
                            yield s, t, n
    else:
        ii = field(-1).sqrt()
        for n in (one, -one, ii, -ii):
            levels = (
                (field(0), field(3) / 2, -field(3) / 2)
                if n * n == one
                else [K for x in (ii, -ii) for K in square_roots(x / 2)]
            )
            for K in levels:
                for N in square_roots(K * K - 1 - n * n):
                    for s, t in lifted_coordinates(K, n * N):
                        yield s, t, n


def tensor_windows(gamma, field):
    def component(i, j, k, u):
        h = i ^ k
        if h != j ^ u:
            return field.zero()
        return sum(
            (
                gamma[h + 4 * b] * (-1) ** (dot(b, k) + dot((0, 3, 2, 1)[b], u))
                for b in range(4)
            ),
            field.zero(),
        )

    windows = []
    for start in range(4):
        side = []
        for h in range(4):
            rows = []
            for k in range(4):
                row = {}
                for j in range(4):
                    word = (j, j ^ h, k, k ^ h)
                    value = component(*(word[-start:] + word[:-start]))
                    if not value.is_zero():
                        row[j] = Scalar(value)
                rows.append(row)
            side.extend(
                (h, {(j, j ^ h): x for j, x in row.items()})
                for row in echelon(rows).values()
            )
        assert len(side) == 6
        windows.append(side)
    return windows


def plane_equations(group):
    K, V = sp.symbols("K V")
    U = K * K
    b, c = (K - 1) / 2, (K + 1) / 2
    assert {(p[8], p[4]) for p, e in group} == {
        (1, 2),
        (1, 3),
        (1, 10),
        (1, 11),
        (8, 4),
        (8, 5),
        (8, 12),
        (8, 13),
        (9, 6),
        (9, 7),
        (9, 14),
        (9, 15),
    }
    for p, e in group:
        assert (e[4] - e[5]) % 2 == 0
    first = (U**2 + 14 * U + 16 * V - 7) ** 2 - 16 * U * (3 * U - 1) ** 2
    second = (U**2 - 2 * U + 16 * V + 9) ** 2 - 16 * U * (U - 3) ** 2
    third = (V - 4 * (U + 1)) ** 2 - 64 * U
    assert (
        sp.expand(first - 256 * (b**4 - 4 * c * b * b + V) * (c**4 + 4 * b * c * c + V))
        == 0
    )
    assert (
        sp.expand(
            second - 256 * (b**4 + 4 * c * b * b + V) * (c**4 - 4 * b * c * c + V)
        )
        == 0
    )
    assert sp.expand(third - (V - 16 * c * c) * (V - 16 * b * b)) == 0


def verify():
    group = template_group()
    plane_equations(group)
    field = fq_default_ctx(17, 16, "u")
    one = field.one()
    ii = field(-1).sqrt()
    phase = [ii**k for k in range(4)]
    for octic in (False, True):
        seen = set()
        reps = 0
        for s, t, n in exceptional_parameters(field, octic):
            if tuple(map(str, (s, t, n))) in seen:
                continue
            gamma = amplitudes(s, t / s, n, one)
            for perm, exp in group:
                row = [gamma[j] * phase[e] for j, e in zip(perm, exp)]
                seen.add(tuple(map(str, recover(row))))
            windows = tensor_windows(gamma, field)
            try:
                assert graded_profile(windows, Scalar(one))
            except AssertionError:
                assert verify_central_profile(windows, Scalar(one))
            reps += 1
        print(
            "Octic" if octic else "Coordinate ends",
            reps,
            "residual orbits certified",
            flush=True,
        )


if __name__ == "__main__":
    verify()
