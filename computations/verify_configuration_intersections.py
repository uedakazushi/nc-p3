"""Restrict the second Veronese ideal to the first and check cyclic chart hypotheses."""

from itertools import combinations, combinations_with_replacement, product

from flint import nmod_mat
import sympy as sp

from verify_klein_elliptic_completion import rational_example
from verify_klein_full_rank_configurations import partner
from verify_klein_two_doublets import modular_points
from verify_unrestricted_doublet_deformations import K, exchange, symmetric_example

PAIRS = list(combinations_with_replacement(K, 2))
QUARTICS = list(combinations_with_replacement(K, 4))


def kernel(matrix, prime):
    padded, size = matrix.nullspace()
    return nmod_mat(
        matrix.ncols(),
        size,
        [int(padded[i, j]) for i in range(matrix.ncols()) for j in range(size)],
        prime,
    )


def squares(points, prime):
    return nmod_mat([[q[i] * q[j] % prime for i, j in PAIRS] for q in points], prime)


def change_of_square_coordinates(a, c, prime):
    first, third = squares(a, prime), squares(c, prime)
    normal = kernel(first.transpose(), prime).transpose()
    equations = nmod_mat(
        [
            [int(normal[i, k]) * int(third[k, j]) % prime for k in range(20)]
            for i, j in product(range(10), repeat=2)
        ],
        prime,
    )
    scales = kernel(equations, prime)
    assert scales.ncols() == 1 and all(scales[i, 0] for i in range(20))
    scaled = nmod_mat(
        [
            [int(scales[i, 0]) * int(third[i, j]) % prime for j in range(10)]
            for i in range(20)
        ],
        prime,
    )
    reduced, _ = first.transpose().rref()
    pivots = [next(i for i in range(20) if reduced[j, i]) for j in range(10)]

    def pivot_rows(matrix):
        return nmod_mat([[int(matrix[i, j]) for j in range(10)] for i in pivots], prime)

    change = pivot_rows(first).inv() * pivot_rows(scaled)
    assert change.det() and first * change == scaled
    return change


def intersection(a, c, prime):
    change = change_of_square_coordinates(a, c, prime)
    quartic_index = {mon: i for i, mon in enumerate(QUARTICS)}

    def multiply(i, j, k, l):
        s = PAIRS.index(tuple(sorted((i, j))))
        t = PAIRS.index(tuple(sorted((k, l))))
        row = [0] * 35
        for u, v in product(range(10), repeat=2):
            index = quartic_index[tuple(sorted(PAIRS[u] + PAIRS[v]))]
            row[index] += int(change[u, s]) * int(change[v, t])
        return row

    restrictions = []
    for i, j in combinations(K, 2):
        for k, l in combinations(K, 2):
            f = multiply(i, k, j, l)
            g = multiply(i, l, j, k)
            restrictions.append([(x - y) % prime for x, y in zip(f, g)])
    return nmod_mat(restrictions, prime).rank()


def rotate(w):
    return sp.Matrix(
        16, 16, lambda ij, kl: w[4 * (kl % 4) + ij // 4, 4 * (ij % 4) + kl // 4]
    )


def rotated_chart(w):
    for mat, expected in [(rotate(w), -sp.Rational(1, 8)), (exchange(rotate(w)), -2)]:
        blocks = [
            sp.Matrix(4, 4, lambda i, k: mat[4 * i + (i ^ h), 4 * k + (k ^ h)])
            for h in K
        ]
        assert [b.rank() for b in blocks] == [2, 2, 1, 1]
        assert all(
            b[i : i + 2, j : j + 2].det()
            for b in blocks[:2]
            for i, j in product((0, 2), repeat=2)
        )
        thin, other = blocks[2:]
        assert all(thin) and all(other)
        ratio = (thin[2, 0] * thin[3, 0] * other[0, 0] * other[1, 0]) / (
            other[2, 0] * other[3, 0] * thin[0, 0] * thin[1, 0]
        )
        assert ratio == expected


if __name__ == "__main__":
    w, _, _ = rational_example()
    rotated_chart(w)
    for mat, prime in [(w, 2657), (symmetric_example(), 12619)]:
        n = exchange(mat)
        a = modular_points(mat, prime)
        b = partner(a, mat, False, prime)
        c = partner(b, n, True, prime)
        d = partner(c, mat.T, False, prime)
        ranks = (intersection(a, c, prime), intersection(b, d, prime))
        assert ranks == ((15, 15) if prime == 2657 else (0, 15))
        print(f"Prime {prime}: restricted Veronese quadrics have ranks {ranks}.")
