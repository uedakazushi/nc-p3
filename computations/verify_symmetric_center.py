"""Reversal on the center, its degree-eight generator and the transverse product defect."""

from functools import cache, reduce
from itertools import combinations, combinations_with_replacement, product
from operator import xor

from flint import nmod_mat
import numpy as np
import sympy as sp

from verify_configuration_intersections import kernel
from verify_quartic_center_order import CentralFiber, integral_quartics
from verify_klein_elliptic_completion import four_relations
from verify_unrestricted_doublet_deformations import (
    K,
    symmetric_example,
    tangent,
    normal_data,
)

PRIME = 1009


def unit_reduce(rows, modulus):
    """Eliminate only unit pivots over Z/p or Z/p², retaining the residual rows."""
    table = rows.copy() % modulus
    pivots = []
    for j in range(table.shape[1]):
        i = len(pivots)
        candidates = np.flatnonzero(table[i:, j] % PRIME)
        if not len(candidates):
            continue
        k = i + candidates[0]
        table[[i, k]] = table[[k, i]]
        table[i] = table[i] * pow(int(table[i, j]), -1, modulus) % modulus
        column = table[:, j].copy()
        column[i] = 0
        table = (table - column[:, None] * table[i]) % modulus
        pivots.append(j)
    return table, pivots


def free_reduction(rows, modulus):
    table, pivots = unit_reduce(rows, modulus)
    assert not table[len(pivots) :].any()
    free = [j for j in range(rows.shape[1]) if j not in pivots]
    quotient = np.zeros((len(free), rows.shape[1]), dtype=np.int64)
    quotient[:, free] = np.eye(len(free), dtype=np.int64)
    quotient[:, pivots] = -table[: len(pivots), free].T % modulus
    return quotient, free


def graded_paths(relations, modulus, bound):
    basis = {(i, 0, h): [()] if h == 0 else [] for i, h in product(K, repeat=2)}
    prefix = {}

    @cache
    def image(i, word):
        if not word:
            return np.ones(1, dtype=np.int64)
        key = (i, len(word), reduce(xor, word), word[0])
        return prefix[key] @ image((i + 1) % 4, word[1:]) % modulus

    for n in range(1, bound + 1):
        for i, h in product(K, repeat=2):
            blocks = [basis[(i + 1) % 4, n - 1, h ^ j] for j in K]
            words = [(j,) + w for j, block in enumerate(blocks) for w in block]
            offset = np.cumsum([0] + list(map(len, blocks)))
            rows = []
            if n >= 2:
                for rule, _ in relations[i]:
                    grade = reduce(xor, next(iter(rule)))
                    for tail in basis[(i + 2) % 4, n - 2, grade ^ h]:
                        row = np.zeros(len(words), dtype=np.int64)
                        for (j, k), c in rule.items():
                            row[offset[j] : offset[j + 1]] += c * image(
                                (i + 1) % 4, (k,) + tail
                            )
                        rows.append(row % modulus)
            table = np.array(rows, dtype=np.int64).reshape(len(rows), len(words))
            proj, free = free_reduction(table, modulus)
            basis[i, n, h] = [words[j] for j in free]
            for j in K:
                prefix[i, n, h, j] = proj[:, offset[j] : offset[j + 1]]
    for i in K:
        for n in range(bound + 1):
            assert (
                sum(len(basis[i, n, h]) for h in K) == (n + 1) * (n + 2) * (n + 3) // 6
            )
    return basis, image


def center_equations(basis, image, degree, h):
    offsets = np.cumsum([0] + [len(basis[i, degree, h]) for i in K])
    rows = []
    for i, j in product(K, repeat=2):
        k = (i + 1) % 4
        left = np.array([image(i, w + (j,)) for w in basis[i, degree, h]]).T
        right = np.array([image(i, (j,) + w) for w in basis[k, degree, h]]).T
        row = np.zeros((len(left), offsets[-1]), dtype=np.int64)
        row[:, offsets[i] : offsets[i + 1]] = left
        row[:, offsets[k] : offsets[k + 1]] = -right
        rows.append(row)
    return np.vstack(rows), offsets


def sections(basis, degree, h, offsets, vector):
    return [
        {
            w: int(c)
            for w, c in zip(basis[i, degree, h], vector[offsets[i] : offsets[i + 1]])
            if c
        }
        for i in K
    ]


def reversal_generator(relations):
    basis, image = graded_paths(relations, PRIME, 9)
    for degree, signs in [
        (4, [(3, 0), (2, 0), (0, 0), (0, 0)]),
        (8, [(8, 1), (6, 0), (0, 0), (0, 0)]),
    ]:
        for h in K:
            equations, offsets = center_equations(basis, image, degree, h)
            space = kernel(nmod_mat(equations.tolist(), PRIME), PRIME)
            reversal = np.zeros((offsets[-1], offsets[-1]), dtype=np.int64)
            for i in K:
                k = (3 - i) % 4
                reversal[offsets[k] : offsets[k + 1], offsets[i] : offsets[i + 1]] = (
                    np.array([image(k, w[::-1]) for w in basis[i, degree, h]]).T
                )
            reversed_space = nmod_mat(reversal.tolist(), PRIME) * space
            anti = kernel(reversed_space + space, PRIME)
            assert (
                space.ncols() - (reversed_space - space).rank(),
                anti.ncols(),
            ) == signs[h]
            if anti.ncols():
                generator = sections(
                    basis, degree, h, offsets, (space * anti).entries()
                )
    return generator


def linear_combination(identity, z, target):
    a, b, c = (m.entries() for m in (identity, z, target))
    i = next(i for i, v in enumerate(a) if v)
    j = next(j for j in range(len(a)) if a[i] * b[j] - a[j] * b[i])
    pair = nmod_mat([[int(a[t]), int(b[t])] for t in (i, j)], PRIME)
    coeff = pair.inv() * nmod_mat(2, 1, [int(c[t]) for t in (i, j)], PRIME)
    u, v = map(int, coeff.entries())
    assert u * identity + v * z == target
    return u, v


def symmetric_fibers(relations, quartics, anti):
    samples = []
    for k, params in enumerate(
        [(0, 0, 0, 0), (1, 2, 3, 4)] + [(1, j, 1 + 2 * j, 3 - j) for j in range(5)]
    ):
        fiber = CentralFiber(relations, quartics, params, PRIME)
        z = fiber.fifth_action(quartics[4])
        y = nmod_mat(64, 64, PRIME)
        for op, c in zip(fiber.loop_operators(anti[0]), anti[0].values()):
            y += c * op
        identity = fiber.prefix(0, 12, ())
        if k == 0:
            assert [z.rank(), y.rank(), (z * y).rank()] == [32, 2, 1]
        if k == 1:
            assert fiber.center_dimension() == 4 and fiber.trace_rank() == 64
        u = linear_combination(identity, z, z * z)
        v = linear_combination(identity, z, y * y)
        if k >= 2:
            samples.append([*u, *v])
        fiber.prefix.cache_clear()
    s, t, z = variables = sp.symbols("s t z")
    inv = nmod_mat([[j**i for i in range(5)] for j in range(5)], PRIME).inv()
    coeff = inv * nmod_mat(samples, PRIME)

    def form(col, deg):
        assert all(not coeff[j, col] for j in range(deg + 1, 5))
        return sum(int(coeff[j, col]) * s ** (deg - j) * t**j for j in range(deg + 1))

    q = z * z - form(1, 1) * z - form(0, 2)
    branch = form(3, 3) * z + form(2, 4)
    assert (
        nmod_mat(
            [list(map(int, row)) for row in sp.hessian(q, variables).tolist()], PRIME
        ).rank()
        == 3
    )
    jac = sp.Matrix([sp.derive_by_array(f, variables).tolist() for f in (q, branch)])
    equations = [q, branch] + [
        jac[:, [i, j]].det() for i, j in combinations(range(3), 2)
    ]
    for j, x in enumerate(variables):
        chart = {**dict.fromkeys(variables[:j], 0), x: 1}
        eq = [f.subs(chart) for f in equations]
        if j == 2:
            assert any(int(f) % PRIME for f in eq)
        else:
            assert list(sp.groebner(eq, *variables[j + 1 :], modulus=PRIME)) == [1]
    print(
        "Center ranks and signs verified. The branch quartic meets the smooth conic transversely."
    )


def base_quadric(relations, quartics):
    # Interpolate z² = A(x)+B(x)z in all four parameters.
    grid = [
        (1, i, j, k) for i in range(3) for j in range(3 - i) for k in range(3 - i - j)
    ]
    values = []
    pairs = list(combinations_with_replacement(K, 2))
    for point in grid:
        fiber = CentralFiber(relations, quartics, point, PRIME)
        z = fiber.fifth_action(quartics[4])
        values.append(linear_combination(fiber.prefix(0, 12, ()), z, z * z))
        fiber.prefix.cache_clear()
    quad = nmod_mat([[p[i] * p[j] for i, j in pairs] for p in grid], PRIME).inv()
    quad *= nmod_mat([[a] for a, b in values], PRIME)
    pivots = list(sp.Matrix(grid).T.rref()[1])
    lin = nmod_mat([list(grid[i]) for i in pivots], PRIME).inv()
    lin *= nmod_mat([[values[i][1]] for i in pivots], PRIME)
    assert nmod_mat(grid, PRIME) * lin == nmod_mat([[b] for a, b in values], PRIME)
    x = sp.symbols("x0:5")
    form = x[4] ** 2 - sum(int(lin[i, 0]) * x[i] * x[4] for i in K)
    form -= sum(int(quad[k, 0]) * x[i] * x[j] for k, (i, j) in enumerate(pairs))
    assert (
        nmod_mat(
            [list(map(int, row)) for row in sp.hessian(form, x).tolist()], PRIME
        ).rank()
        == 5
    )
    print("The base quadric in P4 has rank five.")


def transverse_product(w):
    modulus = PRIME**2
    project, _ = normal_data(w)
    _, [direction], _ = tangent(w, 0, -1, project)
    flip = [4 * j + i for i, j in product(K, repeat=2)]
    relations = []
    for k, (base, delta) in enumerate(
        [(w, direction), (w, -direction), (w.T, direction.T), (w.T, -direction.T)]
    ):
        indices = list(base.to_DM().rref()[1])
        table = (base + PRIME * delta)[:, indices].T
        if k % 2:
            table = table[:, flip]
        relations.append(
            [
                (
                    {
                        divmod(j, 4): int(c.p) * pow(int(c.q), -1, modulus) % modulus
                        for j, c in enumerate(row)
                        if c
                    },
                    None,
                )
                for row in table.tolist()
            ]
        )
    basis, image = graded_paths(relations, modulus, 8)
    quartics = []
    for h in (0, 1):
        equations, offsets = center_equations(basis, image, 4, h)
        vectors, _ = free_reduction(equations, modulus)
        assert len(vectors) == [3, 2][h]
        quartics.extend((h, sections(basis, 4, h, offsets, v)) for v in vectors)
    for h in (0, 1):
        cols = []
        for (a, u), (b, v) in combinations_with_replacement(quartics, 2):
            if a ^ b == h:
                cols.append(
                    np.concatenate(
                        [
                            sum(
                                (c * d % modulus) * image(i, s + t)
                                for s, c in u[i].items()
                                for t, d in v[i].items()
                            )
                            % modulus
                            for i in K
                        ]
                    )
                )
        reduced, pivots = unit_reduce(np.array(cols).T, modulus)
        residual = reduced[len(pivots) :]
        assert len(pivots) == [8, 6][h] and not (residual % PRIME).any()
        assert bool(residual.any()) == (h == 0)
    print(
        "Quartic products have rank 14 on the divisor and a nonzero transverse derivative."
    )


if __name__ == "__main__":
    w = symmetric_example()
    rules, quartics = integral_quartics(four_relations(w, w), PRIME)
    symmetric_fibers(rules, quartics, reversal_generator(rules))
    base_quadric(rules, quartics)
    transverse_product(w)
