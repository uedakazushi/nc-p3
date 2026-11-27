"""Two-plane blocks in a nontrivially filtered rational change of a volume tensor."""
from itertools import permutations, product
from math import prod

from sympy import Matrix, eye

pairs = list(product(range(4), repeat=2))
shapes = list(product(range(2), repeat=4))


def sign(p):
    return (-1) ** sum(a > b for i, a in enumerate(p) for b in p[i+1:])


# The first two basis vectors are normal. In the original volume every term
# has two normals; the four independent shears add higher-weight corrections.
changes = []
for i in range(4):
    g = eye(4)
    g[0, 1] = i + 1
    g[2, 3] = -2 * i
    g[0, 2] = 3 * i + 1
    g[1, 3] = 1 - i
    changes.append(g)
w = {p: sum(sign(q) * prod(changes[i][p[i], q[i]] for i in range(4))
            for q in permutations(range(4))) for p in product(range(4), repeat=4)}
count = lambda p: sum(j < 2 for j in p)
assert {count(p) for p, coeff in w.items() if coeff} == {2, 3, 4}
initial = {p: value if count(p) == 2 else 0 for p, value in w.items()}
for i in range(4):
    def at(p):
        return initial[p[-i:] + p[:-i]] if i else initial[p]
    flat = Matrix([[at(row + col) for col in pairs] for row in pairs])
    mixed = [[pair for pair in pairs if tuple(int(j < 2) for j in pair) == pattern]
             for pattern in ((1, 0), (0, 1))]
    assert flat.rank() == 6
    for rows, cols in product(mixed, repeat=2):
        assert Matrix([[at(row + col) for col in cols] for row in rows]).det() != 0
    block = {bits: at((bits[0], bits[1], bits[2]+2, bits[3]+2)) for bits in shapes}
    invariant = sum((-1)**sum(p) * block[p] * block[tuple(1-j for j in p)] for p in shapes)
    assert abs(invariant) == 4
print("PASS: six blocks, both complementary lines, all mixed determinants and four Levi invariants")
