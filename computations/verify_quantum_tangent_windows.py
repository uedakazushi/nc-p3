"""Check all 512 supports, two fourth-root bases, and the AS3 linear parameter computation."""

from itertools import product

import sympy as sp

from quantum_tangent_windows import (
    equivalent,
    finite_edges,
    local_lengths,
    period,
    permuted,
)
from ternary_ore_parameters import parameters


def quantum(p, q, r):
    basis = sp.eye(9)
    return sp.Matrix.hstack(
        *(
            basis[:, 3 * a + b] - v * basis[:, 3 * b + a]
            for a, b, v in [(0, 1, p), (1, 2, q), (2, 0, r)]
        )
    )


def main():
    for triple, dim in [
        ((2, 3, 5), 0),
        ((sp.Rational(1, 8), 2, 2), 6),
        ((1, 1, 2), 2),
        ((-1, -1, -1), 9),
        ((sp.I, sp.I, sp.I), 9),
    ]:
        _, _, f = parameters(quantum(*triple))
        assert f.cols == dim
        if dim == 9:
            zeta = triple[0]
            for root in {zeta, zeta**2, zeta**3}:
                expected = sum(zeta**j == root for j in (1, 2, 3)) * 3
                assert 9 - (f - root * sp.eye(9)).rank() == expected
    example = ((2, 3, 5), (7, 11, 13), (17, 19, 23))
    for zeta, n in [(-1, 2), (sp.I, 4)]:
        assert local_lengths(zeta, example) == (14, (1,) * 6)
        assert period(zeta, example) == n
        target, d = permuted(sp.sympify(zeta), example, (1, 0, 2))
        assert equivalent(zeta, example, target, d)
        gauge = (sp.Integer(2), sp.Integer(3), sp.Integer(5))
        moved = [
            [
                v * gauge[(a + 1) % 3] ** s * gauge[(a + 2) % 3] ** (2 - s) / gauge[a]
                for s, v in enumerate(row)
            ]
            for a, row in enumerate(example)
        ]
        assert equivalent(zeta, example, zeta, moved)
    assert not equivalent(-1, example, sp.I, example)
    # Discriminant zero at one torus point, simple roots on the two remaining lines.
    assert local_lengths(-1, ((1, 2, -1), example[1], example[2])) == (
        14,
        (1, 1, 1, 1, 2),
    )
    assert period(sp.I, ((2, 0, 5), (7, 0, 13), (17, 0, 23))) == 2
    assert period(sp.I, ((1, 0, 0),) * 3) == 1
    counts = {}
    for flat in product(range(2), repeat=9):
        rows = [flat[3 * i : 3 * i + 3] for i in range(3)]
        if finite_edges(rows):
            length, _ = local_lengths(sp.I, rows)
            counts[length] = counts.get(length, 0) + 1
    assert counts == dict(zip(range(14, 21), (8, 24, 54, 50, 36, 6, 2)))
    print("finite support counts by local length", sorted(counts.items()))
    print(
        "PASS quantum resonance dimensions, shift eigenspaces, periods, relabellings, degenerations"
    )


if __name__ == "__main__":
    main()
