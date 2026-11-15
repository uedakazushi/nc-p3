#!/usr/bin/env python3
"""Exact rational relation-insertion ranks; finite-degree regression, not AS proof."""
from fractions import Fraction as F
from itertools import product
from sympy import Matrix, eye, kronecker_product


def insertion_columns(windows, start, degree):
    for pos in range(degree - 1):
        for left in product(range(4), repeat=pos):
            for right in product(range(4), repeat=degree - pos - 2):
                for relation in windows[(start + pos) % 4]:
                    col = {}
                    for pair, value in relation.items():
                        index = 0
                        for letter in left + pair + right:
                            index = 4 * index + letter
                        col[index] = F(value)
                    yield col


def exact_rank(columns):
    """Sparse Gaussian elimination over Q, with no numerical tolerance."""
    pivots = {}
    for column in columns:
        col = dict(column)
        while col:
            pivot = min(col)
            if pivot not in pivots:
                scale = col[pivot]
                pivots[pivot] = {i: x / scale for i, x in col.items()}
                break
            scale = col[pivot]
            for i, x in pivots[pivot].items():
                value = col.get(i, F(0)) - scale * x
                if value:
                    col[i] = value
                else:
                    col.pop(i, None)
    return len(pivots)


def seven_window(sign):
    # Order W,Y,Z,X; u=1,v=2,r=3,s=4,a=1,b=2,c=3.
    return [
        {(1, 0): 1, (0, 1): 1},
        {(2, 0): 1, (0, 2): 1},
        {(2, 1): 1, (1, 2): 1, (0, 0): -1},
        {(3, 0): 1, (0, 3): 1, (1, 1): -1, (2, 2): -2,
         (1, 2): -sign, (0, 0): F(sign, 2)},
        {(3, 1): 1, (1, 3): 1, (0, 0): -1, (2, 2): -1,
         (1, 1): -3, (0, 1): -sign, (0, 2): -2 * sign},
        {(3, 2): 1, (2, 3): 1, (1, 1): -4, (0, 1): -3 * sign},
    ]


def dense_columns(columns, rows):
    cols = list(columns)
    matrix = Matrix.zeros(rows, len(cols))
    for j, col in enumerate(cols):
        for i, value in col.items():
            matrix[i, j] = value
    return matrix


def relation_matrix(relations):
    return Matrix(16, 6, lambda i, j: relations[j].get(divmod(i, 4), 0))


def quintic_d3_rank(windows, start):
    # W3 = (R_i tensor V) intersection (V tensor R_{i+1}).
    first = kronecker_product(relation_matrix(windows[start]), eye(4))
    second = kronecker_product(eye(4), relation_matrix(windows[(start + 1) % 4]))
    syzygies = first.row_join(-second).nullspace()
    assert len(syzygies) == 4
    coordinates = Matrix.hstack(*(v[:24, :] for v in syzygies))
    # A quotient map V^2 -> C_2 and a section, without choosing normal words.
    tail = relation_matrix(windows[(start + 3) % 4])
    quotient2 = Matrix.hstack(*tail.T.nullspace()).T
    pivots = quotient2.rref()[1]
    assert len(pivots) == 10
    section = Matrix.zeros(16, 10)
    inverse = quotient2[:, list(pivots)].inv()
    for row, pivot in enumerate(pivots):
        section[pivot, :] = inverse[row, :]
    assert quotient2 * section == eye(10)
    # Quotient of the last three letters by both relation positions.
    relations3 = dense_columns(insertion_columns(windows, (start + 2) % 4, 3), 64)
    quotient3 = Matrix.hstack(*relations3.T.nullspace()).T
    assert quotient3.rows == 20
    differential = (kronecker_product(eye(6), quotient3)
                    * kronecker_product(coordinates, section))
    assert differential.shape == (120, 40)
    return exact_rank(({i: F(differential[i, j]) for i in range(120)
                        if differential[i, j]} for j in range(40)))


def main():
    commutative = [{(b, a): 1, (a, b): -1}
                   for a in range(4) for b in range(a + 1, 4)]
    samples = {'polynomial': [commutative] * 4,
               'seven_signed': [seven_window((-1) ** i) for i in range(4)]}
    for name, windows in samples.items():
        for start in range(4):
            ranks = [exact_rank(insertion_columns(windows, start, n))
                     for n in (2, 3, 4, 5)]
            assert ranks == [6, 44, 221, 968], (name, start, ranks)
            print(f'PASS {name} start={start}: ranks={ranks}, h=(10,20,35,56)')
        for start in range(4):
            rank = quintic_d3_rank(windows, start)
            assert rank == 36, (name, start, rank)
            print(f'PASS {name} start={start}: quintic d3 rank={rank}')
        dual = []
        for relations in windows:
            rows = [[r.get((a, b), 0) for a in range(4) for b in range(4)]
                    for r in relations]
            dual.append([{divmod(j, 4): F(x) for j, x in enumerate(v) if x}
                         for v in Matrix(rows).nullspace()])
        for start in range(4):
            dims = [4 ** n - exact_rank(insertion_columns(dual, start, n))
                    for n in (2, 3, 4)]
            assert dims == [6, 4, 1], (name, start, dims)
            print(f'PASS {name} start={start}: dual dimensions={dims}')
    # Negative control outside the Frobenius locus: six forbidden pairs
    # among the first three letters; it must fail the rank-221 filter.
    bad = [{(a, b): 1} for a in range(3) for b in range(a, 3)]
    rank = exact_rank(insertion_columns([bad] * 4, 0, 4))
    assert rank != 221
    print(f'PASS negative control: rank={rank}; no Frobenius claim for this input')


if __name__ == '__main__':
    main()
