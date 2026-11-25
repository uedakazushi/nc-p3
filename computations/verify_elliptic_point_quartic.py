"""Exact quartic multiplication on two diagonal quadrics in a Sklyanin example."""
from itertools import product
from math import comb

from flint import fmpq, fmpq_mat


def words(length):
    return product(range(4), repeat=length)


def word_index(word):
    value = 0
    for letter in word:
        value = 4 * value + letter
    return value


def ideal_rows(relations, degree):
    rows = []
    for offset in range(degree - 1):
        for left in words(offset):
            for right in words(degree - 2 - offset):
                for relation in relations.values():
                    row = [0] * (4 ** degree)
                    for pair, coeff in relation.items():
                        row[word_index(left + pair + right)] = coeff
                    rows.append(row)
    return rows


def central_diagonals(relations):
    rows = ideal_rows(relations, 3)
    reduced, rank = fmpq_mat(rows).rref()
    pivots = [next(j for j in range(64) if reduced[i, j]) for i in range(rank)]
    columns = []
    for a in range(4):
        values = []
        for b in range(4):
            row = [0] * 64
            row[word_index((a, a, b))] += 1
            row[word_index((b, a, a))] -= 1
            for i, pivot in enumerate(pivots):
                if row[pivot]:
                    scalar = row[pivot]
                    row = [x - scalar * reduced[i, j] for j, x in enumerate(row)]
            values.extend(row)
        columns.append(values)
    reduced, rank = fmpq_mat(columns).transpose().rref()
    assert rank == 2
    pivots = [next(j for j in range(4) if reduced[i, j]) for i in range(rank)]
    basis = []
    for free in sorted(set(range(4)) - set(pivots)):
        vector = [0] * 4
        vector[free] = 1
        for i, pivot in enumerate(pivots):
            vector[pivot] = -reduced[i, free]
        basis.append(vector)
    return basis


def main():
    relations = {}
    for (i, j, k), alpha in zip(((1, 2, 3), (2, 3, 1), (3, 1, 2)),
                                (2, 3, fmpq(-5, 7))):
        relations[i, 0] = {(0, i): 1, (i, 0): -1, (j, k): -alpha, (k, j): -alpha}
        relations[i, 1] = {(0, i): 1, (i, 0): 1, (j, k): -1, (k, j): 1}
    quadrics = central_diagonals(relations)
    rows = ideal_rows(relations, 4)
    rank = fmpq_mat(rows).rank()
    assert 256 - rank == 35
    products = []
    for left, right in product(quadrics, repeat=2):
        row = [0] * 256
        for a, b in words(2):
            row[word_index((a, a, b, b))] = left[a] * right[b]
        products.append(row)
    assert fmpq_mat(rows + products).rank() == rank + 3
    skew = [a - b for a, b in zip(products[1], products[2])]
    assert fmpq_mat(rows + [skew]).rank() == rank
    def section_dim(d):
        return 4 * d if d > 0 else 1

    for n in range(50):
        assert sum((h + 1) * section_dim(n - 2 * h)
                   for h in range(n // 2 + 1)) == comb(n + 3, 3)
    print("PASS: two central quadrics; unique quartic relation is alternating and invertible")


if __name__ == "__main__":
    main()
