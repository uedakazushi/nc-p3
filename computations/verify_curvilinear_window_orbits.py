"""Exact checks of the finite-window orbit coordinates over the algebraic closure of Q."""

import sympy as sp

from curvilinear_window_orbits import (
    equivalent,
    labelled_key,
    reorder_factors,
    torus_invariants,
)
from verify_finite_square_chart import frame


def jet_change(m, coefficients, unit):
    """Coefficient matrix of f(t) -> u(t) f(h(t)) modulo t^m."""
    t = sp.symbols("t")
    h = sum(v * t**i for i, v in enumerate(coefficients))
    u = sum(v * t**i for i, v in enumerate(unit))
    columns = [sp.Poly(u * h**i, t) for i in range(m)]
    return sp.Matrix(m, m, lambda i, j: columns[j].nth(i))


def test_jet():
    lengths = (4,) + (1,) * 16
    # The simple points form a labelled projective frame on every side.
    inputs = tuple(
        sp.eye(4).col_join(
            sp.Matrix([[1, s, s**2, s**3] for s in range(1 + side, 17 + side)])
        )
        for side in range(3)
    )
    changed = []
    for side, b in enumerate(inputs):
        local = jet_change(4, [0, 2, 3, -1], [3 + side, -1, side, 1])
        matrix = sp.diag(local, *[side + i + 2 for i in range(16)])
        basis = sp.eye(4)
        basis[0, 3], basis[2, 1] = 3, 2 + side
        changed.append(matrix * b * basis)
    key = labelled_key(lengths, inputs)
    assert labelled_key(lengths, changed) == key
    independent = sp.diag(jet_change(4, [0, 2], [1]), sp.eye(16))
    assert labelled_key(lengths, (inputs[0], independent * inputs[1], inputs[2])) != key
    moved_lengths, moved = reorder_factors(lengths, changed, list(range(1, 17)) + [0])
    assert equivalent(lengths, inputs, moved_lengths, moved)
    # A small local-first example with nonconsecutive pivots and several moving rows.
    lengths = (6,)
    b = sp.Matrix(
        [
            [1, 2, 0, 1],
            [0, 0, 0, 0],
            [1, 0, 0, 0],
            [1, 2, 0, 1],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
    planes = [b, b[::-1, :], sp.Matrix(6, 4, lambda i, j: (i + 2) ** j)]
    changes = [jet_change(6, [0, 3, 2, 0, -1, 1], [2, -1, s + 1]) for s in range(3)]
    assert labelled_key(lengths, planes) == labelled_key(
        lengths, [c * p for c, p in zip(changes, planes)]
    )


def main():
    x, y = frame()
    lengths = (1,) * 20
    original = (x, y, x)
    changed = []
    for side, b in enumerate(original):
        basis = sp.eye(4)
        basis[0, 1], basis[2, 3] = 3, side + 1
        changed.append(sp.diag(*range(side + 1, side + 21)) * b * basis)
    assert labelled_key(lengths, original) == labelled_key(lengths, changed)
    test_jet()
    weights = [[2, 0], [0, 2], [1, 1]]
    a, b, c = map(sp.Integer, (2, 3, 5))
    key = torus_invariants(weights, [a, b, c])
    assert key == torus_invariants(weights, [4 * a, 9 * b, 6 * c])
    assert key != torus_invariants(weights, [a, b, 2 * c])
    print(
        "PASS: reduced AS example, nonreduced common jets, relabelling, distinct relative jets and integer weight kernel"
    )


if __name__ == "__main__":
    main()
