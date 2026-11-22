"""Explicit diagonal changes of generators at all four ramification disks over 101."""

from flint import nmod
from rank_fourteen_function_field import relations
from sympy import QQ_I
from verify_klein_center_lift import echelon, field_paths
from verify_rank_fourteen_boundary import (
    flattenings,
    homogeneous_tensor,
    ramified_coefficients,
    specialize,
)
from verify_rank_fourteen_geometric_reduction import verify_central_profile

WEIGHTS = {
    1: ((-9, -9, 1, -17), (4, -14, 12, 12), (2, 2, -6, -24), (-15, 3, -7, -7)),
    -1: ((-7, -7, 3, -15), (-3, -21, 5, 5), (-3, -3, 7, -11), (-23, -5, 3, 3)),
    -3: ((10, 10, -7, 11), (-12, 6, -13, -13), (-1, -1, 0, 18), (3, -15, 2, 2)),
    None: ((8, 8, -9, 9), (-12, 6, -13, -13), (9, 9, -8, 10), (13, -5, -6, -6)),
}


def residue(q):
    return nmod(int(q.numerator), 101) / nmod(int(q.denominator), 101)


def unramified_limits(tensor):
    i = QQ_I(0, 1)
    t = -1 + 2 * i
    windows = relations()
    for component in (i, -i):
        for s in (1, -1):
            y = (t + 1) * (t + 3) * s
            values = specialize(tensor, t, component, y)
            assert len(values) == 64
            for z in values.values():
                assert all(residue(z.x) + r * residue(z.y) for r in (10, 91))
            rules = [
                [(h, specialize(row, t, component, y)) for h, row in side]
                for side in windows
            ]
            assert all(
                len(
                    echelon(
                        [
                            {4 * a + b: c for (a, b), c in row.items() if c}
                            for _, row in side
                        ]
                    )
                )
                == 6
                for side in rules
            )
            basis, _ = field_paths(rules, QQ_I.one)
            for p in range(4):
                assert [
                    sum(len(basis[p, n, h]) for h in range(4)) for n in (3, 4, 5)
                ] == [20, 38, 72]


if __name__ == "__main__":
    tensor = homogeneous_tensor()
    unramified_limits(tensor)
    for t0, weights in WEIGHTS.items():
        data = ramified_coefficients(tensor, t0)
        face = {}
        for ix, (order, z, sheet) in data.items():
            excess = 9 * order - sum(weights[p][a] for p, a in enumerate(ix))
            assert excess >= 0
            if excess == 0:
                assert sheet == 0
                face[ix] = z
        assert len(face) == 32
        for i in (10, 91):
            values = {ix: residue(z.x) + i * residue(z.y) for ix, z in face.items()}
            rules, opposite = flattenings(values)
            assert list(map(len, rules)) == [6] * 4
            assert all(
                len({ix[p] for ix, z in values.items() if z}) == 4 for p in range(4)
            )
            assert verify_central_profile(rules, nmod(1, 101))
            print(
                t0, i, "opposite", opposite, "standard four-cut through 13", flush=True
            )
