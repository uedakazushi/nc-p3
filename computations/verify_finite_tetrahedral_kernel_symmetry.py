"""Verify the identities governing the four reduced rank-one directions."""
from itertools import combinations

from sympy import Matrix, eye, symbols, zeros

vectors = list(eye(3).columnspace())
vectors.append(-sum(vectors, zeros(3, 1)))
assert all(Matrix.hstack(*triple).det() != 0 for triple in combinations(vectors, 3))
x, y, z = symbols("x y z")
left = Matrix([x, y, z])
right = Matrix(symbols("r:3"))
normal = left.cross(right)
assert left.dot(normal).expand() == 0
for v in vectors:
    assert (v.cross(normal) - left * v.dot(right) + right * v.dot(left)).applyfunc(
        lambda p: p.expand()).is_zero_matrix
forms = [v.dot(left) for v in vectors]
for n in (3, 4):
    product = 1
    for f in forms[:n]:
        product *= f
    for f in forms[:n]:
        assert (product / f).cancel().is_polynomial(x, y, z)
# Off-diagonal entries of the four-point pencil have two distinct roots.
for i, j in combinations(range(4), 2):
    assert Matrix([[1, i], [1, j]]).det() != 0
print("PASS: exterior product, tetrahedral products and interpolation on the line")
