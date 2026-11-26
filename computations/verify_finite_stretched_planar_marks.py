"""Check the stretched multiplication table and the minor identity over Z."""
from sympy import Matrix, eye, symbols, zeros

for e in range(2, 13):
    # Basis 1, x, ..., x^e, y.
    x, y = zeros(e+2), zeros(e+2)
    for i in range(e):
        x[i+1, i] = 1
    y[e+1, 0] = y[e, e+1] = 1
    assert x*y == y*x == zeros(e+2)
    assert y**2 == x**e and x**(e+1) == zeros(e+2)
    assert x**2 * y == zeros(e+2) and y**3 == zeros(e+2)
a = Matrix(3, 3, symbols("a:9"))
m = Matrix(3, 3, symbols("m:9"))
am = a*m
for i in range(3):
    for j in range(3):
        k = (j+1) % 3
        expansion = sum(a[k, h]*(m[h, k]*m[i, j]-m[h, j]*m[i, k]) for h in range(3))
        assert (expansion-am[k, k]*m[i, j]+am[k, j]*m[i, k]).expand() == 0
print("PASS: stretched multiplication and the two-minor expansion")
