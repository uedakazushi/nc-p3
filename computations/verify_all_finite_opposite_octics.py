"""Check the two nonreduced grids by exact quartic identities, and enumerate their support graphs."""

from itertools import combinations, product

import sympy as sp


def support_graphs():
    for sizes in ((4,), (1, 3), (2, 2), (1, 1, 2), (1, 1, 1, 1)):
        pairs = list(combinations(range(len(sizes)), 2))
        valid = []
        for lengths in product(*(range(sizes[i] * sizes[j] + 1) for i, j in pairs)):
            degrees = [sum(v for pair, v in zip(pairs, lengths) if i in pair)
                       for i in range(len(sizes))]
            if degrees == [2*n for n in sizes]:
                valid.append(lengths)
        if sizes in ((4,), (1, 3)):
            assert valid == []
        elif sizes == (2, 2):
            assert valid == [(4,)]
        elif sizes == (1, 1, 2):
            assert valid == [(0, 2, 2)]
        else:
            assert len(valid) == 3
            assert all(sorted(v) == [0, 0, 1, 1, 1, 1] for v in valid)


def check_identities():
    a, b, c, d = sp.symbols("a b c d")
    metric = sp.Matrix([[0, 0, 2], [0, -1, 0], [2, 0, 0]])
    matrices = [
        sp.Matrix([[-a, a, 0], [-b, b-a, a], [c, c, 0], [d, d+c, c]]),
        sp.Matrix([[0, a, 0], [0, b, a], [a, c, 0], [b, d+a, c]]),
    ]
    for number, (matrix, cubic) in enumerate(zip(matrices, (a*a*c, a**3))):
        gram = matrix * metric * matrix.T

        def minor(rows, cols):
            return gram.extract(rows, cols).det()

        p, q, r, s, t = (0, 1), (0, 2), (0, 3), (1, 2), (2, 3)
        A, B, F, G = minor(p, q), minor(p, s), minor(p, t), minor(r, r)
        if number == 0:
            cert = (-A/4, (A-B)/4, F/8, -(F+G)/8)
        else:
            cert = (-A/2, -minor(p, r)/2, -minor(q, r)/2, (A+F)/2-G/4)
        assert all(sp.expand(x*cubic-y) == 0 for x, y in zip((a, b, c, d), cert))


if __name__ == "__main__":
    support_graphs()
    check_identities()
    print("Verified the length-four support graphs and all eight quartic certificates over Q.")
