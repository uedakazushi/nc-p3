"""Check the exact cubic projection, division by the bracket, and jet identities."""
from sympy import Matrix, symbols

s, t, u = symbols("s t u")
cubic = Matrix([1, t, t**2, t**3])
# The tangent line at t=0 annihilates the last two coordinates.
assert cubic[2:].count(0) == 0
assert [entry / t**2 for entry in cubic[2:]] == [1, t]
assert list(cubic[1:, :].applyfunc(lambda f: f / t)) == [1, t, t**2]
kernel = Matrix([t**2, -s*t, s**2])
assert (Matrix([[s, t, 0], [0, s, t]]) * kernel).is_zero_matrix
multiply = Matrix([[-t, 0, 0], [1, -t, 0], [0, 1, -t], [0, 0, 1]])
assert multiply[1:, :].det() == 1
assert (Matrix([[1, u, u**2, u**3]]) * multiply
        - (u-t) * Matrix([[1, u, u**2]])).applyfunc(lambda f: f.expand()).is_zero_matrix

rows = [Matrix(1, 3, symbols(f"k{i}_0:3")) for i in range(3)]
source = Matrix([1, t, t**2])
word = source * sum((row * t**i for i, row in enumerate(rows)), Matrix.zeros(1, 3))
jets = [word.applyfunc(lambda f: f.expand().coeff(t, j)) for j in range(3)]
for j, jet in enumerate(jets):
    assert jet[j, :] == rows[0]
    assert jet[j + 1:, :].is_zero_matrix
    for i in range(j):
        assert jet[i, :] == rows[j-i]
print("PASS: tangent projection and constant right factor in the order-three jet plane")
