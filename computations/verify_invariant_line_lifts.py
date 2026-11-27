"""Normal-weight ranks and Schur equations on exact split and filtered windows."""
from itertools import combinations, permutations, product
from math import prod

from sympy import Matrix, zeros


WORDS = tuple(product(range(4), repeat=4))
PAIRS = tuple(product(range(4), repeat=2))


def weight(word):
    return sum(letter < 2 for letter in word)


def flatten(w, start):
    def value(a, b):
        q = a + b
        return w.get(q[-start:] + q[:-start] if start else q, 0)
    return Matrix([[value(a, b) for b in PAIRS] for a in PAIRS])


def tangent_matrix(w, n):
    words = [q for q in WORDS if weight(q) == n]
    rows = []
    for start in (0, 1):
        flat = flatten(w, start)
        left = Matrix.hstack(*flat.T.nullspace()).T
        right = Matrix.hstack(*flat.nullspace())
        cols = []
        for q in words:
            perm = q[start:] + q[:start]
            col = left[:, PAIRS.index(perm[:2])] * right[PAIRS.index(perm[2:]), :]
            cols.append(Matrix(list(col)))
        rows.append(Matrix.hstack(*cols))
    full = Matrix.vstack(*rows)
    return full[[i for i, row in enumerate(full.tolist()) if any(row)], :]


def quantum(coeffs):
    factors = dict(zip(combinations(range(4), 2), coeffs))
    return {q: prod(-factors[q[j], q[i]] for i, j in combinations(range(4), 2)
                    if q[i] > q[j]) for q in permutations(range(4))}


def shear(w):
    lifts = [Matrix([[1+i, 2-i], [-2*i, 3]]) for i in range(4)]
    def entry(source, target):
        if source == target:
            return 1
        if source >= 2 and target < 2:
            return lifts[slot][target, source-2]
        return 0
    out = {}
    for q in WORDS:
        total = 0
        for p, value in w.items():
            factor = 1
            for slot in range(4):
                factor *= entry(p[slot], q[slot])
            total += factor * value
        out[q] = total
    return out


def verify():
    for coeffs, quartic_rank in [([1]*6, 13), ([2, 3, 5, 7, 11, 13], 16)]:
        w = quantum(coeffs)
        cubics = tangent_matrix(w, 3)
        quartics = tangent_matrix(w, 4)
        assert 40 <= cubics.rank() <= 48
        assert quartics.rank() == quartic_rank
        lifted = shear(w)
        for start in (0, 1):
            base = flatten(w, start)
            cubic = flatten({q:v for q,v in lifted.items() if weight(q) == 3}, start)
            top = flatten({q:v for q,v in lifted.items() if weight(q) == 4}, start)
            cols, rows = base.rref()[1], base.T.rref()[1]
            inv = base.extract(rows, cols).inv()
            general_inverse = zeros(16)
            for j, b in enumerate(cols):
                for k, a in enumerate(rows):
                    general_inverse[b,a] = inv[j,k]
            left = Matrix.hstack(*base.T.nullspace()).T
            right = Matrix.hstack(*base.nullspace())
            assert left * cubic * right == zeros(10)
            assert left * (top - cubic * general_inverse * cubic) * right == zeros(10)
            assert (base + cubic + top).rank() == 6
        print("cubic constraint rank", cubics.rank(), "quartic constraint rank", quartic_rank)


if __name__ == '__main__':
    verify()
