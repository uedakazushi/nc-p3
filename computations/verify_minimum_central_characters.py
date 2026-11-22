"""Separate the twenty central characters by lifting labeled cycles along a rational tangent."""

from itertools import product
from math import prod
from random import Random

import numpy as np
import sympy as sp

from verify_klein_elliptic_completion import rational_example
from verify_klein_full_rank_configurations import partner
from verify_klein_two_doublets import modular_points
from verify_unrestricted_doublet_deformations import K, exchange, normal_data, tangent
from verify_symmetric_center import (
    graded_paths,
    center_equations,
    free_reduction,
    unit_reduce,
)

PRIME = 2657
MODULUS = PRIME**2


def residue(x):
    return int(x.p) * pow(int(x.q), -1, MODULUS) % MODULUS


def lifted_relations(w, reverse):
    first, _ = normal_data(w)
    second, _ = normal_data(reverse)

    def project(t):
        return first(t).col_join(second(exchange(t)))

    rng = Random(8103)
    direction = sp.zeros(16)
    for h in K:
        _, vectors, _ = tangent(w, h, 0, project)
        for t in vectors:
            direction += rng.randrange(1, PRIME) * t
    assert project(direction).is_zero_matrix
    near = w + PRIME * direction
    other = exchange(near)
    flip = [4 * j + i for i, j in product(K, repeat=2)]
    tables = []
    relations = []
    for i, (base, perturbed) in enumerate(
        ((w, near), (reverse, other), (w.T, near.T), (reverse.T, other.T))
    ):
        cols = list(base.to_DM().rref()[1])
        mat = perturbed[:, cols].T
        if i % 2:
            mat = mat[:, flip]
        table = np.array(
            [[residue(c) for c in row] for row in mat.tolist()], dtype=np.int64
        )
        tables.append(table)
        relations.append(
            [
                ({divmod(j, 4): int(c) for j, c in enumerate(row) if c}, None)
                for row in table
            ]
        )
    return tables, relations


def lift_cycle(point, tables):
    indices = [next(j for j in K if row[j] % PRIME) for row in point]
    point = np.array(
        [
            row * pow(int(row[j]), -1, MODULUS) % MODULUS
            for row, j in zip(point, indices)
        ]
    )
    free = [(i, j) for i, j in product(K, repeat=2) if j != indices[i]]

    def values(coords):
        return np.concatenate(
            [
                tables[i]
                @ (np.outer(coords[i], coords[(i + 1) % 4]).ravel() % MODULUS)
                % MODULUS
                for i in K
            ]
        )

    errors = values(point)
    assert not (errors % PRIME).any()
    columns = []
    for i, j in free:
        shifted = point.copy()
        shifted[i, j] += 1
        columns.append((values(shifted) - errors) % PRIME)
    augmented = np.column_stack([np.array(columns).T, (-errors // PRIME) % PRIME])
    reduced, pivots = unit_reduce(augmented, PRIME, PRIME)
    assert pivots == list(range(12))
    for (i, j), v in zip(free, reduced[:12, -1]):
        point[i, j] += PRIME * v
    assert not values(point).any()
    return point


def main():
    w, _, reverse = rational_example()
    tables, relations = lifted_relations(w, reverse)
    basis, image = graded_paths(relations, MODULUS, 5, prime=PRIME, characters=False)
    equations, offsets = center_equations(basis, image, 4, 0)
    vectors, _ = free_reduction(equations, MODULUS, PRIME)
    assert len(vectors) == 5 and not (equations @ vectors.T % MODULUS).any()
    a = modular_points(w, PRIME)
    b = partner(a, w, False, PRIME)
    c = partner(b, reverse, True, PRIME)
    d = partner(c, w.T, False, PRIME)
    words = basis[0, 4, 0]
    points = []
    for cycle in zip(a, b, c, d):
        coords = lift_cycle(np.array(cycle, dtype=np.int64), tables)
        evaluations = [
            prod(int(coords[i, j]) for i, j in enumerate(word)) % MODULUS
            for word in words
        ]
        value = vectors[:, : offsets[1]] @ np.array(evaluations) % MODULUS
        scale = next(int(v) for v in value if v % PRIME)
        points.append(tuple(int(v) * pow(scale, -1, MODULUS) % MODULUS for v in value))
    assert len({tuple(v % PRIME for v in q) for q in points}) == 8
    assert len(set(points)) == 20
    print(
        "Eight distinct specialized projective characters; all twenty separate modulo 2657²."
    )


if __name__ == "__main__":
    main()
