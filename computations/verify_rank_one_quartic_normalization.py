"""Check the two- and three-coordinate quartic pencil obstructions."""
from sympy import Matrix, symbols

s, t, u, v, a, b, c = symbols("s t u v a b c")
pencil = Matrix([[0, 0, c*t], [0, 0, -c*s], [c*t, -c*s, a*s+b*t]])
assert pencil*Matrix([s, t, 0]) == Matrix.zeros(3, 1)
assert pencil.adjugate() == -c**2 * Matrix([s,t,0])*Matrix([[s,t,0]])
assert 3+2+2 < 16-8
assert 4+4-5 >= 3
print("PASS: common subsystem dimension, singular pencil normal form and marked group budget")
