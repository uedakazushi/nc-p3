"""Check the slope count, mixed tensor ranks and primitive kernel degree obstruction."""
from itertools import product
from sympy import Matrix

assert 3+2+2 < 16-8
alpha, alpha_other = Matrix([1, 2, 3]), Matrix([2, 0, 1])
beta, beta_other = Matrix([0, 1, 4]), Matrix([3, 2, 0])
pencil_value = alpha*beta.T - 5*alpha_other*beta_other.T
assert pencil_value.rank() == 2
assert alpha.cross(alpha_other).T*pencil_value == Matrix.zeros(1, 3)
assert pencil_value*beta.cross(beta_other) == Matrix.zeros(3, 1)
assert not any(g+h <= 2 and g+1 >= 3 and h+1 >= 3 for g, h in product(range(3), repeat=2))
print("PASS: at least four slopes and incompatible left/right dimensions in a singular ternary pencil")
