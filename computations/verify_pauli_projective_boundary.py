"""Classify regularity on the complement of the projective rank-fifteen Pauli chart."""

from concurrent.futures import ProcessPoolExecutor
from itertools import product

import sympy as sp
from flint import fq_default_ctx

from flint_rational_function import rational_function_field
from verify_klein_center_lift import echelon, field_paths
from verify_rank_fifteen_formal_boundary import dual_relations
from verify_rank_fifteen_parameter_orbits import amplitudes, dot, pattern, template_group
from verify_rank_fifteen_pencil_linear_residues import graded_profile
from verify_rank_fifteen_regular_orbits import square_roots, tensor_windows
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile


def relations(gamma, start):
    def component(a, b, c, d):
        if a ^ c != b ^ d:
            return 0
        return sum(
            gamma[(a ^ c) + 4*v] * (-1)**(dot(v, c) + dot((0, 3, 2, 1)[v], d))
            for v in range(4)
        )

    side = []
    for h in range(4):
        rows = []
        for k in range(4):
            row = {}
            for j in range(4):
                w = (j, j ^ h, k, k ^ h)
                value = component(*(w[-start:] + w[:-start]))
                if value:
                    row[j] = value
            rows.append(row)
        side.extend((h, {(j, j ^ h): v for j, v in row.items()})
                    for row in echelon(rows).values())
    return side


def rank_nine(job):
    bits, sign = job
    field = rational_function_field(("u",))
    u = field.from_sympy(sp.Symbol("u"))
    a, c = (u + 1/u)/2, (u - 1/u)/2
    pairs = []
    for root, bit in zip((sign*a*c, a, c), bits):
        pairs.extend((field.zero, root) if bit else (root, field.zero))
    gamma = pattern(field.one, c*c, *pairs)
    sides = [relations(gamma, p) for p in range(3)]
    dimensions = []
    for h in range(4):
        rows = []
        for p, side in enumerate(sides):
            for color, relation in side:
                for j in range(4):
                    other = (j, j ^ h ^ color)
                    row = {}
                    for pair, coefficient in relation.items():
                        word = other[:p] + pair + other[p:]
                        column = sum(a * 4**k for k, a in enumerate(word))
                        row[column] = coefficient
                    rows.append(row)
        dimensions.append(64 - len(echelon(rows)))
    assert dimensions == [10, 10, 8, 8]


def rank_four():
    for sign, swap in product((sp.S.One, -sp.S.One), (False, True)):
        L, M = (0, sign) if swap else (sign, 0)
        gamma = pattern(0, sp.S.One, L, M, 0, 0, 0, 0)
        assert [len(relations(gamma, p)) for p in range(4)] == [4]*4


def exceptional_rank_fourteen():
    field = fq_default_ctx(101, 4, "e")
    one = field.one()
    unit = Scalar(one)
    group = template_group()
    # Some signed permutation makes any specified one of a,b,a+b the new a.
    assert {perm[8] for perm, _ in group} == {1, 8, 9}
    phases = [field(-1).sqrt()**j for j in range(4)]

    def key(row):
        return tuple(str(c / row[8]) for c in row)

    for level, expected in ((1, (256, 2)), (2, (128, 1))):
        seen = set()
        count = representatives = 0
        for K, H in product((one, -one), (field(level), field(-level))):
            rs = [K + H + u for u in square_roots((K + H)**2 - 1)]
            ws = [K - H + u for u in square_roots((K - H)**2 - 1)]
            ns = [n for n2 in square_roots(-H**2) for n in square_roots(n2)]
            for r, w in product(rs, ws):
                for s, t, n in product(square_roots(r), square_roots(w), ns):
                    count += 1
                    row = amplitudes(s, t/s, n, one)
                    if key(row) in seen:
                        continue
                    representatives += 1
                    for perm, exp in group:
                        transformed = [phases[e]*row[j] for j, e in zip(perm, exp)]
                        if not transformed[8].is_zero():
                            seen.add(key(transformed))
                    windows = tensor_windows(row, field)
                    dual, _ = field_paths(dual_relations(windows, unit), unit)
                    for p in range(4):
                        assert [sum(len(dual[p, d, h]) for h in range(4))
                                for d in range(6)] == [1, 4, 6, 4, 1, 0]
                    assert sum(c != field.zero() for c in row) == 14
                    try:
                        assert graded_profile(windows, unit)
                    except AssertionError:
                        assert verify_central_profile(windows, unit)
        assert (count, representatives) == expected


if __name__ == "__main__":
    rank_four()
    with ProcessPoolExecutor(max_workers=4) as pool:
        list(pool.map(rank_nine, product(product(range(2), repeat=3), (1, -1))))
    exceptional_rank_fourteen()
    print("All 16 rank-nine charts fail in degree four; rank-four and rank-fourteen boundaries verified.")
