"""Square cover, projective diagonal equivalence, two-shift and the quadratic dual."""

from itertools import permutations, product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_decomp


def phase(i):
    return (-1) ** ((i & 1) * (i >> 1))


def verify():
    labels = range(4)
    words = [(i, i ^ b, i ^ b ^ e, i ^ e) for i, b, e in product(labels, range(3), range(3))]
    colors = list(product(range(3), repeat=2))
    weights = sp.Matrix([[int(w[p] == j) for p, j in product(labels, repeat=2)] for w in words])
    pairs = sp.Matrix([[int((w[0] ^ w[1], w[1] ^ w[2]) == pair) for pair in colors] for w in words])
    smith, moves, _ = smith_normal_decomp(weights, domain=sp.ZZ)
    rank = weights.rank()
    assert [smith[j, j] for j in range(rank)] == [1]*11 + [2, 2]
    tests = moves[rank:, :] * pairs
    assert {tuple(int(z % 2) for z in row) for row in tests.tolist()} == {
        (0,) * 9, tuple(int(b != e) for b, e in colors)
    }
    for i, j, k, t in words:
        assert phase(i) * phase(j) * phase(k) * phase(t) == phase(i ^ k)

    for perm in permutations((1, 2, 3)):
        linear = (0, *perm)
        assert all(linear[x ^ y] == linear[x] ^ linear[y] for x, y in product(labels, repeat=2))
        def beta(x, y):
            return ((x >> 1)*(y & 1) + (linear[x] >> 1)*(linear[y] & 1)) % 2
        def refinement(x):
            return (beta(1, 1)*(x & 1) + beta(2, 2)*(x >> 1) + 2*beta(1, 2)*(x & 1)*(x >> 1)) % 4
        for x, y in product(labels, repeat=2):
            assert (refinement(x ^ y) - refinement(x) - refinement(y) - 2*beta(x, y)) % 4 == 0
        assert sum(refinement(b ^ e) for b, e in colors if b != e) % 4 == 0
        for d in labels:
            assert (refinement(d) + d.bit_count() + linear[d].bit_count()) % 2 == 0
    assert sum((b >> 1)*((b ^ e) & 1) for b, e in colors if b != e) % 2 == 1

    z = sp.symbols("z:4")
    block = sp.Matrix(4, 4, lambda i, k: (-1)**((i >> 1)*((i ^ k) & 1))*z[i ^ k])
    metric = sp.diag(1, -1, -1, 1)
    norm = (sp.Matrix(z).T * metric * sp.Matrix(z))[0]
    assert sp.simplify(block.T * metric * block - norm * metric) == sp.zeros(4)
    assert sp.simplify(block * metric * block.T - norm * metric) == sp.zeros(4)


if __name__ == "__main__":
    verify()
    print("PASS: signed parameter quotient, simultaneous permutations, two-shift and dual cut")
