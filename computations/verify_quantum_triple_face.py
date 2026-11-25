"""Exact checks of the nine-parameter cyclic tensor and its attraction cones."""

from itertools import combinations, product

import sympy as sp


def outer(left, right):
    return {i + j: a * b for i, a in left.items() for j, b in right.items()}


def basis_tensors():
    q = sp.symbols("Q")
    background = {}
    corrections = []
    for a, b, c in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
        commutator = {(0, a): 1, (a, 0): -1}
        terms = ({(b, c): 1, (c, b): -q}, {(a, a): 1},
                 {(a, 0): sp.S.Half, (0, a): sp.S.Half}, {(0, 0): 1})
        for s, term in enumerate(terms):
            tensor = {}
            for left, right in ((commutator, term), (term, commutator)):
                for word, value in outer(left, right).items():
                    tensor[word] = tensor.get(word, 0) + value
            if s == 0:
                background.update(tensor)
            else:
                corrections.append({w: v for w, v in tensor.items() if v})
    return background, corrections


def weight(word):
    return tuple(word.count(i) - 1 for i in range(4))


def verify():
    background, corrections = basis_tensors()
    assert len(background) == 24
    for tensor in [background, *corrections]:
        assert all(tensor.get(w[1:] + w[:1]) == -v for w, v in tensor.items())
    weights = []
    for tensor in corrections:
        chars = {weight(w) for w in tensor}
        assert len(chars) == 1
        weights.append(chars.pop())
    assert len(set(weights)) == 9
    assert [sum(weights[i][j] for i in (0, 3, 6)) for j in range(4)] == [0] * 4
    for missing in range(3):
        cocharacter = [6, 1, 1, 1]
        cocharacter[missing + 1] = -2
        degrees = [sum(a * b for a, b in zip(w, cocharacter)) for w in weights]
        assert [i for i, d in enumerate(degrees) if d <= 0] == [3 * missing]
        for a in range(3):
            expected = [1, 0, 0, 0]
            expected[a + 1] = -1
            assert [x - y for x, y in zip(weights[3*a + 2], weights[3*a + 1])] == expected

    # Verify the chart determinant at arbitrary symbolic coefficients by its triangular blocks.
    pairs = list(product(range(4), repeat=2))
    cyclic = [(1, 2, 3), (2, 3, 1), (3, 1, 2)]
    rows = [(0, a) for a, _, _ in cyclic] + [(b, c) for _, b, c in cyclic]
    cols = rows[3:] + rows[:3]
    scalars = sp.symbols("z:9")
    for rotation in range(2):
        matrix = sp.zeros(16)
        for scale, tensor in zip([1, *scalars], [background, *corrections]):
            for word, v in tensor.items():
                w = word[rotation:] + word[:rotation]
                matrix[pairs.index(w[:2]), pairs.index(w[2:])] += scale * v
        minor = matrix.extract([pairs.index(i) for i in rows], [pairs.index(j) for j in cols])
        assert minor.det() == 1

    # Fourier squares generate the exceptional ternary relation space over Q(omega).
    omega = sp.Symbol("w")
    for i in range(3):
        for a, b, c in cyclic:
            assert sp.rem(omega**(2*i*a)-omega**(i*(b+c)), omega**2+omega+1, omega) == 0
    for p, q in combinations(range(3), 2):
        assert any((p - e) % 3 == (q + e) % 3 for e in range(3))


if __name__ == "__main__":
    verify()
    print("PASS: all nine cyclic perturbations, two unit chart minors, three support planes and exceptional Fourier relations")
