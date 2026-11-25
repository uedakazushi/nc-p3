"""Labelled orbit coordinates for three four-planes in a curvilinear point algebra.

Work in characteristic zero. Integer powers in torus invariants test equivalence
over the algebraic closure, even when the displayed coordinates are rational.
"""

from itertools import accumulate
from math import factorial

import sympy as sp


def chart(matrix):
    """Use the earliest independent rows as pivots for the column space."""
    matrix = sp.Matrix(matrix)
    pivots = matrix.T.rref()[1]
    if len(pivots) != 4:
        raise ValueError("each coordinate space must be four-dimensional")
    return pivots, matrix * matrix[list(pivots), :].inv()


def jet_generators(lengths):
    """A generator consists of three actions on the underlying point algebra."""
    size = sum(lengths)
    generators = []
    for start, m in zip(accumulate((0, *lengths[:-1])), lengths, strict=True):
        for k in range(1, m):
            shift = sp.SparseMatrix(
                size, size, {(start + j + k, start + j): 1 for j in range(m - k)}
            )
            for side in range(3):
                generators.append(
                    tuple(
                        shift if j == side else sp.SparseMatrix.zeros(size)
                        for j in range(3)
                    )
                )
            if k < m - 1:
                derivative = sp.SparseMatrix(
                    size, size, {(start + j + k, start + j): j for j in range(1, m - k)}
                )
                generators.append((derivative,) * 3)
    return generators


def act(generator, scalar, planes, pivots):
    result = []
    for x, plane, rows in zip(generator, planes, pivots, strict=True):
        term = plane
        moved = plane
        for i in range(1, plane.rows):
            term = x * term
            if term.is_zero_matrix:
                break
            moved = moved + scalar**i / factorial(i) * term
        result.append(moved * moved[list(rows), :].inv())
    return result


def unipotent_chart(lengths, inputs):
    if len(inputs) != 3 or any(a.shape != (sum(lengths), 4) for a in inputs):
        raise ValueError(
            "expected three matrices with four columns and one row per basis vector"
        )
    if not lengths or min(lengths) < 1:
        raise ValueError("local lengths must be positive")
    pivots, planes = map(list, zip(*(chart(a) for a in inputs)))
    generators = jet_generators(lengths)
    for side, rows in enumerate(pivots):
        for r in range(sum(lengths)):
            if r in rows:
                continue
            for j in reversed(range(4)):
                if rows[j] >= r or not generators:
                    continue
                b = planes[side]

                def velocity(g):
                    x = g[side]
                    return (x.row(r) * b[:, j])[0] - sum(
                        b[r, i] * (x.row(p) * b[:, j])[0]
                        for i, p in enumerate(rows)
                        if b[r, i]
                    )

                speeds = [velocity(g) for g in generators]
                lead = next((i for i, s in enumerate(speeds) if s), None)
                if lead is None:
                    continue
                move = generators[lead]
                if b[r, j]:
                    planes = act(move, -b[r, j] / speeds[lead], planes, pivots)
                    assert planes[side][r, j] == 0
                # Speeds depend only on earlier coordinates, which the move fixed.
                generators = [
                    tuple(
                        a - s / speeds[lead] * x for a, x in zip(g, move, strict=True)
                    )
                    if s
                    else g
                    for i, (g, s) in enumerate(zip(generators, speeds, strict=True))
                    if i != lead
                ]
    return pivots, planes


def torus_data(lengths, pivots, planes):
    positions = [(a, d) for a, m in enumerate(lengths) for d in range(m)]
    support, weights, values = [], [], []
    for side, (rows, b) in enumerate(zip(pivots, planes, strict=True)):
        for r, (a, d) in enumerate(positions):
            if r in rows:
                continue
            for j, p in enumerate(rows):
                if b[r, j] == 0:
                    continue
                other, e = positions[p]
                weight = [0] * (4 * len(lengths))
                weight[4 * a + side] += 1
                weight[4 * other + side] -= 1
                weight[4 * a + 3] += d
                weight[4 * other + 3] -= e
                support.append((side, r, j))
                weights.append(weight)
                values.append(b[r, j])
    return tuple(support), weights, values


def torus_invariants(weights, values):
    """Remove each torus parameter by invertible integer row operations."""
    weights = [row[:] for row in weights]
    values = list(values)
    for j in range(len(weights[0]) if weights else 0):
        pivot = next((i for i, row in enumerate(weights) if row[j]), None)
        if pivot is None:
            continue
        for i, row in enumerate(weights):
            if i == pivot or row[j] == 0:
                continue
            w = weights[pivot]
            a, b = w[j], row[j]
            u, v, d = map(int, sp.gcdex(a, b))
            h, k = a // d, b // d
            weights[pivot] = [u * x + v * y for x, y in zip(w, row, strict=True)]
            weights[i] = [h * y - k * x for x, y in zip(w, row, strict=True)]
            left, right = values[pivot], values[i]
            values[pivot], values[i] = left**u * right**v, left ** (-k) * right**h
        del weights[pivot], values[pivot]
    return tuple(values)


def labelled_key(lengths, planes):
    pivots, normalized = unipotent_chart(lengths, planes)
    support, weights, values = torus_data(lengths, pivots, normalized)
    return tuple(lengths), tuple(pivots), support, torus_invariants(weights, values)


def reorder_factors(lengths, planes, order):
    """Transport a tuple along a permutation, retaining the local jet bases."""
    if sorted(order) != list(range(len(lengths))):
        raise ValueError("expected a permutation of the local factors")
    starts = list(accumulate((0, *lengths[:-1])))
    rows = [starts[a] + j for a in order for j in range(lengths[a])]
    return tuple(lengths[a] for a in order), tuple(b[rows, :] for b in planes)


def equivalent(source_lengths, source, target_lengths, target):
    """Try length-compatible identifications lazily, with prefix-rank pruning."""
    if sorted(source_lengths) != sorted(target_lengths):
        return False
    source_key = labelled_key(source_lengths, source)

    def signatures(lengths, planes):
        starts = list(accumulate((0, *lengths[:-1])))
        return [
            tuple(
                b[start : start + j, :].rank() for b in planes for j in range(1, m + 1)
            )
            for start, m in zip(starts, lengths, strict=True)
        ]

    src = signatures(source_lengths, source)
    dst = signatures(target_lengths, target)

    def match(order, unused):
        if not unused:
            return (
                labelled_key(*reorder_factors(target_lengths, target, order))
                == source_key
            )
        a = len(order)
        return any(
            match(order + (b,), unused - {b})
            for b in sorted(unused)
            if source_lengths[a] == target_lengths[b] and src[a] == dst[b]
        )

    return match((), set(range(len(target_lengths))))
