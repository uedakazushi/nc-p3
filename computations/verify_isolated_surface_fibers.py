"""Recover the linear inverse on two projected quadrics and check the Euler identity."""
from sympy import Matrix, Poly, expand, symbols

x, y, z = symbols("x y z")
parameters = symbols("a:12")
linear_inverse = Matrix(3, 4, parameters)
# Projections from a smooth point on a smooth quadric and on a cone.
for coordinates in ([x*x, x*y, x*z, y*z], [x*x, x*y, y*y, x*z]):
    image = linear_inverse * Matrix(coordinates)
    minors = Matrix([x, y, z]).cross(image)
    equations = [coeff for f in minors for coeff in Poly(f, x, y, z).coeffs()]
    solutions = Matrix(equations).jacobian(parameters).nullspace()
    assert len(solutions) == 1
    inverse = solutions[0].reshape(3, 4)
    assert inverse.rank() == 3
    assert all(expand(f) == 0 for f in Matrix([x,y,z]).cross(inverse * Matrix(coordinates)))

n = symbols("n")
weights = (1, -4, 6, -4, 1)
assert expand(sum(c * (n-i+2) * (n-i+1) / 2 for i, c in enumerate(weights))) == 0
print("PASS: unique linear inverse in both quadric models, and zero plane Euler polynomial")
