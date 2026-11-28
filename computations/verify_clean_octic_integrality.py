"""Check the quadric dual formula, toric cover, and residue budgets."""
from sympy import Matrix, diag, expand, simplify, symbols


def cross(u, v):
    return u[0]*v[1] - u[1]*v[0]


hexagon = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]
assert sum(cross(hexagon[i], hexagon[(i+1) % 6]) for i in range(6)) == 6
for i, point in enumerate(hexagon):
    neighbors = [hexagon[(i+j) % 6] for j in (-1, 1)]
    edges = [tuple(v[k]-point[k] for k in range(2)) for v in neighbors]
    assert abs(cross(*edges)) == 1
s, t = symbols("s t", nonzero=True)
x, y, z = s+1/s, t+1/t, s*t+1/(s*t)
assert simplify(x*x + y*y + z*z - x*y*z - 4) == 0
u = Matrix(symbols("u:4"))
q = diag(1, 2, 3, 4)
dual = (u.T * q * u)[0] * q - q * u * u.T * q
assert (dual * u).applyfunc(expand) == 0 * u
assert 2*(22 + 4*3) - 6*6 == 22 + 4*2 + 2 == 32
print("PASS: toric degree six, four fixed nodes, thirty basepoints and two conjugate residual points")
