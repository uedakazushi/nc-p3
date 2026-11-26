"""Exact identities in the pencil and binary-form arguments of the proof."""
from itertools import combinations

from sympy import Matrix, eye, kronecker_product, symbols, zeros

s, t, u, v = symbols("s t u v")
for d in (1, 2, 3):
    syzygies = zeros(d + 1, d)
    for j in range(d):
        syzygies[j, j] = -t
        syzygies[j + 1, j] = s
    forms = Matrix([[u ** (d - j) * v ** j for j in range(d + 1)]])
    expected = Matrix([[(s * v - t * u) * u ** (d - 1 - j) * v ** j
                        for j in range(d)]])
    assert (forms * syzygies - expected).applyfunc(lambda f: f.expand()).is_zero_matrix
    minors = [syzygies.extract(rows, range(d)).det()
              for rows in combinations(range(d + 1), d)]
    assert set(minor.as_coeff_Mul()[1] for minor in minors) == {
        s ** (d - j) * t ** j for j in range(d + 1)
    }

bracket = Matrix([[0, 1], [-1, 0]])
assert kronecker_product(bracket, eye(3)).rank() == 6
for a in (2, 3, 4):
    for b in (2, 3, 4):
        for r in range(min(3, a, b)):
            assert 6 + a + b - r >= 9
print("PASS: restricted pencil identities, two-line adjacent isomorphism and rank bound")
