"""Check the discriminant budgets, the tetrahedral graph and the edge determinant."""
from itertools import combinations
from sympy import Matrix, expand, symbols

assert max(3, 4, 4, 5) < 6
vertices = range(4)
edges = list(combinations(vertices, 2))
for size in range(4, 7):
    for selected in combinations(edges, size):
        degrees = [sum(v in e for e in selected) for v in vertices]
        if max(degrees) <= 2:
            assert size == 4 and degrees == [2]*4
            assert len(set(edges) - set(selected)) == 2
a, b, l, c, d = symbols("a b l c d")
m = Matrix([[a*a + l*c, a*b], [a*b, b*b + l*d]])
assert expand(m.det() - l*(d*a*a+c*b*b+l*c*d)) == 0
assert 10+6+3 < 20 and 10+3+2+4 < 20 and 2*5+1+4 < 16
print("PASS: discriminant counts, all four-vertex graphs, edge determinant and terminal marked bounds")
