"""Flagged finite-window invariants in products of truncated plane or line algebras.

Each input factor is (embedding dimension, smallest n with maximal_ideal**n=0).
For plane factors y spans the chosen cotangent line. The ideal matrix is a column
basis of the full kernel; each coordinate matrix lifts four independent classes.
"""

from itertools import accumulate

import sympy as sp

from curvilinear_window_orbits import torus_invariants, triangular_chart


def monomials(dimension, n):
    if dimension == 0 and n == 1:
        return [(0, 0)]
    if dimension == 1 and n >= 2:
        return [(j, 0) for j in range(n)]
    if dimension == 2 and n >= 3:
        return [(j - b, b) for j in range(n) for b in range(j + 1)]
    raise ValueError("invalid truncated factor")


def first_relations(ideal, mons):
    """Degree-d components of the ideal elements of order at least d."""
    x, y = sp.symbols("x y")
    for d in range(2, sum(mons[-1]) + 1):
        lo = next(i for i, ex in enumerate(mons) if sum(ex) == d)
        same = [i for i, ex in enumerate(mons) if sum(ex) == d]
        kernel = ideal[:lo, :].nullspace()
        if not kernel:
            continue
        leading = (ideal[same, :] * sp.Matrix.hstack(*kernel)).columnspace()
        if leading:
            forms = [
                sum(c * x ** mons[i][0] * y ** mons[i][1] for c, i in zip(vec, same))
                for vec in leading
            ]
            if len(forms) == 1:
                return forms[0]
            if len(forms) == 2:
                f, g = forms
                return sp.expand(f.diff(x) * g.diff(y) - f.diff(y) * g.diff(x))
            raise ValueError("expected a planar complete intersection")
    raise ValueError("no visible initial relation space")


def flagged_key(types, ideal, coordinates):
    blocks = [monomials(*spec) for spec in types]
    n = sum(map(len, blocks))
    if (
        ideal.rows != n
        or any(a.shape != (n, 4) for a in coordinates)
        or len(coordinates) != 3
    ):
        raise ValueError("incompatible ambient dimensions")
    inputs = [ideal] + [ideal.row_join(a) for a in coordinates]
    zero = sp.SparseMatrix.zeros(n)
    generators = []
    weights = [[] for _ in inputs]
    for block, start, mons in zip(
        range(len(blocks)),
        accumulate([0] + [len(b) for b in blocks[:-1]]),
        blocks,
        strict=True,
    ):
        indices = {ex: i + start for i, ex in enumerate(mons)}
        for a, b in mons:
            for side in range(4):
                w = [0] * (5 * len(blocks))
                w[5 * block : 5 * block + 2] = [a, b]
                if side:
                    w[5 * block + 1 + side] = 1
                weights[side].append(w)
        for a, b in mons[1:]:
            multiplication = {}
            for (i, j), col in indices.items():
                dest = (i + a, j + b)
                if dest in indices:
                    multiplication[indices[dest], col] = 1
            mult = sp.SparseMatrix(n, n, multiplication)
            for side in range(1, 4):
                generators.append(tuple(mult if s == side else zero for s in range(4)))
            for var in range(types[block][0]):
                if a + b == 1 and (a, b, var) != (0, 1, 0):
                    continue
                derivative = {}
                for (i, j), col in indices.items():
                    dest = (i + a - int(var == 0), j + b - int(var == 1))
                    value = (i, j)[var]
                    if value and dest in indices:
                        derivative[indices[dest], col] = value
                if derivative:
                    generators.append((sp.SparseMatrix(n, n, derivative),) * 4)
    pivots, normalized = triangular_chart(inputs, generators)
    support, chars, values = [], [], []
    for s, (rows, b) in enumerate(zip(pivots, normalized)):
        for r in range(n):
            if r in rows:
                continue
            for j, p in enumerate(rows):
                if b[r, j]:
                    support.append((s, r, j))
                    chars.append([a - c for a, c in zip(weights[s][r], weights[s][p])])
                    values.append(b[r, j])
    return tuple(types), tuple(pivots), tuple(support), torus_invariants(chars, values)
