"""Adjacent and opposite flattenings in a Walsh basis, including the rank-14 boundary."""

from itertools import product

import sympy as sp
from flint import fq_default_ctx

from verify_klein_center_lift import field_paths
from verify_rank_fifteen_formal_boundary import dual_relations
from verify_rank_fifteen_parameter_orbits import amplitudes, pattern, recover, template_group
from verify_rank_fifteen_pencil_linear_residues import graded_profile
from verify_rank_fifteen_regular_orbits import lifted_coordinates, square_roots, tensor_windows
from verify_rank_fourteen_geometric_reduction import Scalar


def sign(index, character):
    return (-1) ** (int(index & character).bit_count())


def verify_blocks():
    labels = range(4)
    z = sp.symbols("z:4")
    shifts = [
        sp.Matrix(4, 4, lambda i, k: int(k == i ^ d) * (-1) ** ((i >> 1) * (d & 1)))
        for d in labels
    ]
    assert shifts[1] * shifts[2] == -shifts[2] * shifts[1]
    assert shifts[1] ** 2 == shifts[2] ** 2 == sp.eye(4)
    block = sum((v * a for v, a in zip(z, shifts)), sp.zeros(4))
    norm = sum(sign(i, 3) * z[i] ** 2 for i in labels)
    assert sp.expand(block.det() - norm**2) == 0
    for d, character in product(labels, repeat=2):
        vector = sp.Matrix([sign(i, character) for i in labels])
        circulant = sp.Matrix(4, 4, lambda i, j: z[i ^ j])
        assert circulant * vector == vector * sum(z[b] * sign(b, character) for b in labels)
        for j, e in product(labels, repeat=2):
            b = d ^ e
            assert ((j ^ b) >> 1) * (d & 1) % 2 == (((j >> 1) + ((d ^ e) >> 1)) * (d & 1)) % 2

    a, b, p, m, q, n, r, l = sp.symbols("a b p m q n r l")
    matrix = sp.Matrix([[a + b, p, n, 0], [m, -b, l, 0], [-q, r, -a, 0], [0, 0, 0, 0]])
    gamma = pattern(a, b, (p + m) / 2, (p - m) / 2, (n - q) / 2, (n + q) / 2,
                    (l - r) / 2, (l + r) / 2)
    for i, j, k, t in product(labels, repeat=4):
        d = i ^ k
        if d != j ^ t:
            continue
        walsh_value = sum(
            (-1) ** (((i ^ j) & 1) * s + (u >> 1) * (k >> 1)
                     + ((0, 3, 2, 1)[u] >> 1) * (t >> 1))
            * gamma[(d & 2) + s + 4 * u]
            for s, u in product(range(2), range(d & 1, 4, 2))
        ) / 2
        assert sp.expand(walsh_value - (-1) ** ((i >> 1) * (d & 1)) * matrix[i ^ j, j ^ k]) == 0
    signed = sp.Matrix(3, 3, lambda i, j: sign(i ^ j, 3) * matrix[i, j] ** 2)
    sums = list(signed * sp.ones(3, 1)) + list(signed.T * sp.ones(3, 1))
    equations = [r*r + l*l - 2*a*b, p*p - b*b - r*r, m*m - b*b - l*l,
                 q*q - a*a - r*r, n*n - a*a - l*l]
    domain = (a, b, p, m, q, n, r, l)
    assert sp.groebner(sums, *domain) == sp.groebner(equations, *domain)


def verify_rank_fourteen():
    group = template_group()
    assert {perm[1] for perm, _ in group} == {1, 8, 9}
    field = fq_default_ctx(97, 4, "v")
    one = field.one()
    unit = Scalar(one)
    fourth = field(-1).sqrt()
    for level, expected_count in ((1, 128), (2, 64)):
        points = set()
        for h in (field(level), field(-level)):
            for s, t in lifted_coordinates(one, h):
                for n2 in square_roots(-h*h):
                    for n in square_roots(n2):
                        points.add(tuple(field(v) if isinstance(v, int) else v for v in amplitudes(s, t / s, n, one)))
        assert len(points) == expected_count
        while points:
            gamma = next(iter(points))
            orbit = set()
            ordinary = False
            for perm, exponents in group:
                if perm[1] != 1:
                    continue
                row = [gamma[j] * fourth**e for j, e in zip(perm, exponents)]
                row = tuple(v / row[8] for v in row)
                orbit.add(row)
                s, t, n = recover(row)
                h = (s*s + 1/(s*s) - t*t - 1/(t*t)) / 4
                assert h**2 == field(level)**2
                if not ((1+h) * h * (h+2) * (2-h)).is_zero():
                    ordinary = True
            assert len(orbit) == 64 and orbit <= points
            if level == 1:
                assert ordinary
            else:
                windows = tensor_windows(gamma, field)
                basis, _ = field_paths(dual_relations(windows, unit), unit)
                for start in range(4):
                    assert [sum(len(basis[start, n, h]) for h in range(4)) for n in range(6)] == [1, 4, 6, 4, 1, 0]
                assert graded_profile(windows, unit)
            points -= orbit


if __name__ == "__main__":
    verify_blocks()
    verify_rank_fourteen()
    print("PASS: three flattenings, five quadrics, and both rank-fourteen exceptional levels")
