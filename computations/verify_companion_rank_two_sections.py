"""Frobenius cycles of three finite hyperplane sections of a rank-two curve."""

from functools import cache
from itertools import combinations, permutations, product
from math import gcd
from random import Random

from flint import nmod_mat
import sympy as sp

from verify_klein_elliptic_completion import gram_twist


PRIME = 1009


def monomials(n, degree):
    if n == 1:
        return [(degree,)]
    return [(i, *m) for i in range(degree + 1) for m in monomials(n - 1, degree - i)]


def multiply(first, second):
    result = {}
    for a, u in first.items():
        for b, v in second.items():
            key = tuple(x + y for x, y in zip(a, b))
            result[key] = (result.get(key, 0) + u * v) % PRIME
    return result


def minor(matrix, rows, cols):
    result = {}
    for p in permutations(range(3)):
        sign = (-1) ** sum(p[i] > p[j] for i, j in combinations(range(3), 2))
        a, b, c = (matrix[rows[i]][cols[p[i]]] for i in range(3))
        for key, val in multiply(a, multiply(b, c)).items():
            result[key] = (result.get(key, 0) + sign * val) % PRIME
    return result


def section(relations, seed):
    exps = sum((monomials(4, n) for n in range(6)), [])
    index = {e: i for i, e in enumerate(exps)}
    rng = Random(seed)
    coeff = [[int(i == j) for i in range(5)] for j in range(5)] + [
        [rng.randrange(PRIME) for _ in range(5)]
    ]
    matrix = [
        [
            {
                e: sum(coeff[k][n] * relations[k][4 * i + j] for k in range(6)) % PRIME
                for n, e in enumerate(exps[:5])
            }
            for j in range(4)
        ]
        for i in range(4)
    ]
    triples = list(combinations(range(4), 3))
    cubics = [minor(matrix, row, col) for row, col in product(triples, repeat=2)]
    # Deg-5 leading parts span every form at infinity; eliminate in descending total degree.
    rows = []
    for f, m in product(cubics, exps[:15]):
        row = [0] * 126
        for e, c in f.items():
            key = tuple(x + y for x, y in zip(e, m))
            row[125 - index[key]] = c
        rows.append(row)
    matrix, rank = nmod_mat(rows, PRIME).rref()
    assert rank == 106
    pivots = {next(j for j in range(126) if matrix[i, j]): i for i in range(rank)}
    free = [j for j in range(126) if j not in pivots]
    assert all(sum(exps[125 - j]) <= 4 for j in free)
    multiplication = []
    for j in free:
        e = exps[125 - j]
        col = 125 - index[(e[0] + 1, *e[1:])]
        if col in pivots:
            multiplication.append([-int(matrix[pivots[col], k]) % PRIME for k in free])
        else:
            multiplication.append([int(k == col) for k in free])
    factors = nmod_mat(multiplication, PRIME).charpoly().factor()[1]
    assert all(n == 1 for f, n in factors)
    return tuple(sorted(f.degree() for f, n in factors))


@cache
def orbit_types(degrees):
    if not degrees:
        return {()}
    result = set()
    # Group closed point orbits by the Frobenius orbit of their geometric component.
    for mask in range(1, 1 << len(degrees), 2):
        group = [d for j, d in enumerate(degrees) if (mask >> j) & 1]
        rest = tuple(d for j, d in enumerate(degrees) if not (mask >> j) & 1)
        for count in range(1, gcd(*group) + 1):
            if all(d % count == 0 for d in group):
                block = (count, sum(group) // count)
                result.update(tuple(sorted(t + (block,))) for t in orbit_types(rest))
    return result


def check():
    A = sp.Matrix([[5, 2], [3, 7]])
    B = sp.Matrix([[sp.Rational(1, 6), 2], [3, -sp.Rational(1, 4)]])
    w, _, _ = gram_twist(A, B, -sp.Integer(2), sp.S.One)
    rels = [
        [int(sp.numer(c)) * pow(int(sp.denom(c)), -1, PRIME) % PRIME for c in v]
        for v in w.columnspace()
    ]
    assert nmod_mat(rels, PRIME).rank() == 6
    types = [section(rels, seed) for seed in (3, 9, 11)]
    assert types == [(1, 3, 16), (2, 2, 4, 4, 8), (1, 4, 15)]
    assert set.intersection(*(orbit_types(t) for t in types)) == {
        ((1, 20),),
        ((1, 4), (1, 16)),
    }
    print(
        "Three reduced sections have cycle types",
        types,
        "and leave only degrees 4 + 16 once the diagonal quartic is included.",
    )


if __name__ == "__main__":
    check()
