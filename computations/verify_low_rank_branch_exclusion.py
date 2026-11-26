"""Check the binary kernel, plane projection, and integer corridor constraints."""
from sympy import Matrix, symbols

s, t = symbols("s t")
l2 = Matrix([[s, t, 0], [0, s, t]])
assert l2 * Matrix([t*t, -s*t, s*s]) == Matrix.zeros(2, 1)
# A quadratic matrix with a primitive linear right kernel and rank two.
a = Matrix([[s*t, -s*s, t*t], [t*t, -s*t, s*s], [s*t, -s*s, s*t]])
assert a * Matrix([s, t, 0]) == Matrix.zeros(3, 1)
assert a.subs({s: 2, t: 3}).rank() == 2
# At three moving positions with adjacent centers the second forward defect
# must be positive whenever it is nonnegative.
assert [3*d-4 for d in (1, 2) if 3*d-4 >= 0] == [2]
assert all(6-d > 0 for d in (1, 2))
print("PASS: binary Kronecker conic, linear matrix kernel, and both degree contradictions")
