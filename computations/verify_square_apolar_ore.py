"""Exact coefficient and sample checks for the double-plane Ore normal form."""
from itertools import combinations_with_replacement, permutations, product

from sympy import Matrix, LeviCivita, expand, eye, symbols, zeros


indices = range(3)
triples = list(product(indices, repeat=3))
symmetric_triples = list(combinations_with_replacement(indices, 3))
coefficients = dict(zip(symmetric_triples, symbols("c:10")))
polar = [Matrix(3, 3, lambda i, j: coefficients[tuple(sorted((a, i, j)))])
         for a in indices]
p = Matrix(3, 3, lambda i, j: symbols("p%d%d" % tuple(sorted((i, j)))))
alternating = [Matrix(3, 3, lambda i, j: sum(LeviCivita(i, j, b)*p[b, a]
                                             for b in indices)) for a in indices]


def flat(h, c, i, k, j, l):
    return sum(h[a][i, j]*c[a][k, l] + c[a][i, j]*h[a][k, l]
               - h[a][i, l]*c[a][k, j] - c[a][i, l]*h[a][k, j]
               for a in indices)


# The matrix of complementary minors on wedge-square is a cross product.
pairs = [(1, 2), (2, 0), (0, 1)]
contraction = [sum(p[i, j]*polar[a][i, j] for i, j in product(indices, repeat=2))
               for a in indices]
for u, (i, k) in enumerate(pairs):
    for v, (j, l) in enumerate(pairs):
        expected = sum(LeviCivita(u, v, a)*contraction[a] for a in indices)
        assert (flat(alternating, polar, i, k, j, l) - expected).expand() == 0


def compare_jacobian():
    basis = []
    for i, j in combinations_with_replacement(indices, 2):
        entry = zeros(3)
        entry[i, j] = entry[j, i] = 1
        basis.append(entry)
    # Trilinearity reduces the alternating residual identity to these basis tests.
    for left, middle, right in product(basis, repeat=3):
        actual = sum(LeviCivita(i, j, b)*left[b, a]*middle[i, k]*right[j, l]
                     * polar[a][k, l] for b, i, j in permutations(indices)
                     for a, k, l in triples
                     if left[b, a] and middle[i, k] and right[j, l])
        swapped = sum(LeviCivita(b, i, j)*middle[b, a]*left[i, k]*right[j, l]
                      * polar[a][k, l] for b, i, j in permutations(indices)
                      for a, k, l in triples
                      if middle[b, a] and left[i, k] and right[j, l])
        assert expand(actual + swapped) == 0


def gram_example(diagonal, fail_ci=False):
    c = [zeros(3) for _ in indices]
    if fail_ci:
        for a in indices:
            c[a][a, a] = 1
    else:
        for a, i, j in permutations(indices):
            c[a][i, j] = 1
    h = [Matrix(3, 3, lambda i, j: LeviCivita(i, j, a)*diagonal[a]
                + int(a == i)*(j+1) + int(a == j)*(i+1)) for a in indices]
    values = {}
    for i, j in product(range(4), repeat=2):
        if i*j:
            values[i, j] = Matrix([m[i-1, j-1] for m in h+c])
        else:
            values[i, j] = Matrix([int(i+j == a+1) for a in indices] + [0]*3)
    adj = Matrix.hstack(*values.values())
    assert adj.rank() == 6
    g = zeros(6)
    g[:3, 3:] = g[3:, :3] = eye(3)

    def pairing(i, j, k, l):
        return (values[i, j].T*g*values[k, l])[0]

    for i, j, k, l in product(range(4), repeat=4):
        assert pairing(i, j, k, l) == pairing(i, l, k, j)
    for i, k, l, a in product(indices, repeat=4):
        overlap = sum(c[b][i, a]*h[b][l, k]-h[b][i, k]*c[b][l, a]
                      for b in indices)
        relations = sum(c[b][i, k]*h[b][l, a]-h[b][i, a]*c[b][l, k]
                        for b in indices)
        assert overlap == relations
    opp = list(combinations_with_replacement(range(4), 2))
    w = Matrix([[pairing(i, j, k, l) for j, l in opp] for i, k in opp])
    assert w[0, :] == zeros(1, 10) and w[:, 0] == zeros(10, 1)
    assert w.rank() == (8 if any(diagonal) else 6)
    assert w.nullspace() == w.T.nullspace()
    if fail_ci:
        assert c[0].rank() == 1
    return w.rank()


if __name__ == "__main__":
    compare_jacobian()
    print("Generic signed complementary minors and alternating socle form verified.")
    print("Opposite ranks:", [gram_example(d) for d in [(0, 0, 0), (1, 2, 3),
                                                        (1, 1, 0), (1, 0, 0)]])
    print("Non-CI polar containing a pure relation:", gram_example((0, 0, 0), True))
