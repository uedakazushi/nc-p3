#!/usr/bin/env python3
"""Exact certificates for rigidity of the twenty points and the generic cubic obstruction."""
from itertools import combinations, combinations_with_replacement
from random import Random

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from verify_finite_square_chart import frame, products, tangent_dimension


def verify_plane_certificate():
    x, _ = frame()
    q = products(x, x, symmetric=True)
    rows = set(range(x.rows))
    coordinates = [sp.zeros(4, 1), *sp.eye(4).columnspace()]
    normals = [a-b for a, b in combinations(coordinates, 2)]
    planes = [{i for i, value in enumerate(x*a) if value == 0} for a in normals]
    pairs = list(combinations(range(len(planes)), 2))
    for i, plane in enumerate(planes):
        restriction = q[sorted(plane), :]
        dependency, = restriction.T.nullspace()
        assert len(plane) == 7 and restriction.rank() == 6 and all(dependency)
        assert sp.Matrix.hstack(*(n for j, n in enumerate(normals) if j != i)).rank() == 4
        if i:
            overlap = sorted(rows - planes[0] - plane)
            assert any(all(x[list(four), :].det() for four in combinations(five, 4))
                       for five in combinations(overlap, 5))
    assert not set.intersection(*planes)
    for i, j in pairs:
        assert q[sorted(planes[i] | planes[j]), :].rank() == 9
    triangles = []
    for i, j, k in combinations(range(len(planes)), 3):
        complement = sorted(rows - planes[i] - planes[j] - planes[k])
        if x[complement, :].rank() == 4:
            triangles.append({(i, j), (i, k), (j, k)})
    assert len(triangles) == 80
    for edge in pairs:
        reached = {edge}
        while True:
            enlarged = reached.union(*(t for t in triangles if t & reached))
            if enlarged == reached:
                break
            reached = enlarged
        assert len(reached) == len(pairs)
    assert tangent_dimension(x, q) == 0


def verify_flip_obstructions():
    x, y = frame()
    rel = sp.Matrix.hstack(*products(x, y).nullspace())
    insertion = sp.kronecker_product(rel, sp.eye(4))
    swap = [16*(i % 4) + 4*((i // 4) % 4) + i // 16 for i in range(64)]
    reverse = insertion[swap, :]
    assert (insertion - reverse).rank() == 20
    assert (insertion + reverse).rank() == 24

    tensors = []
    minor_indices = list(combinations(range(4), 2))
    pairs = list(combinations_with_replacement(range(6), 2))
    for i, j in pairs:
        tensor = rel[:, i] * rel[:, j].T
        if i != j:
            tensor += tensor.T
        tensors.append(sp.Matrix([tensor[4*a+b, 4*c+d] - tensor[4*a+d, 4*c+b]
                                  for a, c in minor_indices for b, d in minor_indices]))
    generator, = sp.Matrix.hstack(*tensors).nullspace()
    coefficients = sp.zeros(6)
    for (i, j), value in zip(pairs, generator):
        coefficients[i, j] = coefficients[j, i] = value
    assert coefficients.det() != 0
    return rel, coefficients


def verify_generic_obstruction():
    random = Random(17)
    relations = sp.Matrix(16, 6, lambda i, j: random.randrange(-2, 3))
    insertion = sp.kronecker_product(relations, sp.eye(4))
    swap = [16*(i % 4) + 4*((i // 4) % 4) + i // 16 for i in range(64)]
    reverse = insertion[swap, :]
    assert DomainMatrix.from_Matrix(insertion.row_join(reverse)).convert_to(sp.GF(101)).rank() == 48


if __name__ == "__main__":
    verify_plane_certificate()
    verify_generic_obstruction()
    verify_flip_obstructions()
    print("PASS ten plane circuits, 80 spanning triples, closure from all 45 pairs,")
    print("projective frames, reduced square-root fiber, generic cubic rank 48 modulo 101, and both flip obstructions")
