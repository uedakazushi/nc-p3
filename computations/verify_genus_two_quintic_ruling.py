"""Verify the scroll classes and cyclic coefficient vanishings over the rationals."""
from itertools import product

from sympy import SparseMatrix


def equations(edge):
    forbidden = [(3, 0), (2, 3)]
    swaps = [((3, j + 1), (j, 3)) for j in range(2)]
    if edge % 2:
        forbidden = [pair[::-1] for pair in forbidden]
        swaps = [(a[::-1], b[::-1]) for a, b in swaps]
    other_coefficient = -1 if edge == 3 else 1
    return [[(pair, 1)] for pair in forbidden] + [
        [(a, 1), (b, other_coefficient)] for a, b in swaps
    ]


def main():
    assert [m for m in range(6) if 5 >= 2*m and -m*m+5*m-4 == 2] == [2]
    assert [(a, b) for a in range(1, 5) for b in range(a, 5)
            if a+b == 5 and (a-1)*(b-1) == 2] == [(2, 3)]
    words = list(product(range(4), repeat=4))
    word_index = {word: i for i, word in enumerate(words)}
    entries = {}
    row = 0
    for edge in range(4):
        for u, v in product(range(4), repeat=2):
            for eq in equations(edge):
                for pair, coefficient in eq:
                    cyclic = (*pair, u, v)
                    word = cyclic[-edge:] + cyclic[:-edge] if edge else cyclic
                    entries[row, word_index[word]] = coefficient
                row += 1
    relations = SparseMatrix(row, len(words), entries)
    forced = [i for i, w in enumerate(words) if w.count(3) % 2 or (w[0], w[2]) == (3, 2)]
    unit_rows = SparseMatrix(len(forced), len(words), {(j, i): 1 for j, i in enumerate(forced)})
    assert len(forced) == 126
    assert relations.to_DM().rank() == relations.col_join(unit_rows).to_DM().rank()
    print("PASS: both scroll classes and all 126 forced zero coefficients")


if __name__ == "__main__":
    main()
