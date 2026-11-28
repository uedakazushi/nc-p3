"""Check marked degree bounds and the local normal form on a rank-two plane."""
from math import comb
from sympy import Matrix, expand, symbols

for marks in range(13, 21):
    assert max(37, marks + 4 + 15, 2*marks + 4 + 7) < 3*marks
assert comb(5, 3) < 13
x, y, z, l, q = symbols("x y z l q")
a = Matrix([l, x, y, z])
v = Matrix([1, 0, 0, 0])
m = a*a.T + q*v*v.T
for i in range(1, 4):
    for j in range(1, 4):
        assert expand(m[0, 0]*m[i, j] - m[0, j]*m[i, 0] - q*a[i]*a[j]) == 0
print("PASS: residual quintic partitions, ten-point restriction and rank-one conic ideal")
