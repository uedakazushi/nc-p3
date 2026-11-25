"""Exact checks of plane-graph Ore corrections and singular fixed lines."""

import sympy as sp

from ternary_ore_parameters import parameters


def graph_relations(g):
    rows = []
    inverse = g.inv()
    for i, j in [(0, 1), (1, 2), (2, 0)]:
        rel = sp.zeros(3)
        rel[i, :] = inverse[j, :]
        rel[j, :] = -inverse[i, :]
        rows.append(rel.reshape(9, 1))
    return sp.Matrix.hstack(*rows)


def main():
    cases = [
        (sp.eye(3), 15, 1),
        (sp.diag(1, 1, 3), 3, 3),
        (sp.Matrix([[1, 0, 1], [0, 1, 0], [0, 0, 1]]), 3, 1),
        (sp.diag(1, 2, 3), 0, 1),
        (sp.Matrix([[1, 1, 0], [0, 1, 0], [0, 0, 2]]), 0, 1),
        (sp.Matrix([[1, 1, 0], [0, 1, 1], [0, 0, 1]]), 0, 1),
    ]
    for g, dim, eigenvalue in cases:
        _, _, shift = parameters(graph_relations(g))
        assert shift == eigenvalue * sp.eye(dim)
    print("PASS: recurrence and shift on six plane graph Jordan types")
    x, y, z, u, a, b, c, s = sp.symbols("x y z u a b c s")
    p = sp.Matrix([x, y, z])
    e = sp.Matrix([1, 0, 0])
    linear = sp.Matrix([[s, 0, 1]])
    g = sp.eye(3) + e * linear
    f = a * x**2 + b * x * y + c * y**2
    # PBW compatibility is a polynomial identity for all pairs of generators.
    sigma = g.inv() * p
    delta = g.inv() * e * f
    assert all(sp.cancel(h) == 0 for h in (p - sigma).cross(delta))
    # At the marked point u=1, polarize f in the right variable r.
    row = sp.Matrix([[a * x + b * y / 2, b * x / 2 + c * y, 0]])
    matrix = g - e * row
    rhs = matrix.adjugate() * p
    q = (linear * p)[0] - f
    assert all(sp.expand(h) == 0 for h in p.cross(rhs) + q * p.cross(e))
    print("PASS: symbolic overlap and marked chart, including transvections")


if __name__ == "__main__":
    main()
