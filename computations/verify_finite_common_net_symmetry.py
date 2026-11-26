"""Verify the permutation argument over Q and the four-form block rank."""
from itertools import permutations, combinations

from sympy import Matrix, eye, zeros
from sympy.combinatorics import Permutation


def main():
    # Multilinearization on four formal letters. A row is one tensor value.
    words = list(permutations(range(4)))
    basis = eye(24)

    def value(word):
        return basis.row(words.index(tuple(word)))

    def pair(first, second, rest, slot):
        a, b = rest
        positions = ((first, second, a, b), (a, first, second, b),
                     (a, b, first, second), (first, b, a, second))
        return value(positions[slot])

    equations = []
    for i, j in combinations(range(4), 2):
        for rest in permutations(set(range(4)) - {i, j}):
            doubles = [pair(i, j, rest, k) + pair(j, i, rest, k) for k in range(4)]
            equations.extend(row - doubles[0] for row in doubles[1:])
    constraints = Matrix.vstack(*equations)
    solutions = Matrix.hstack(*constraints.nullspace())
    signs = Matrix([Permutation(list(w)).signature() for w in words])
    assert constraints * signs == zeros(constraints.rows, 1)
    for i, j in ((0, 2), (1, 3)):
        rows = []
        for w in words:
            v = list(w)
            v[i], v[j] = v[j], v[i]
            rows.append(value(w) - value(v))
        assert Matrix.vstack(*rows) * solutions == signs * signs.T * solutions / 12
    pairs = list(combinations(range(4), 2))
    wedge = Matrix([[Permutation(list(a + b)).signature() if set(a).isdisjoint(b)
                     else 0 for b in pairs] for a in pairs])
    assert wedge.rank() == 6
    print("PASS: polarized central-square identities, both opposite swaps, wedge rank six")


if __name__ == "__main__":
    main()
