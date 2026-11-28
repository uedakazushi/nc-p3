"""Symbolic identities of the pencil-kernel cubic and its separating minors."""
from sympy import Matrix, symbols, factor

u,v,b1,b2,b3,a,c,e,g = symbols("u v b1 b2 b3 a c e g")
N = Matrix([[a,u*b1,u*b2,u*b3],
            [u*b1,v*b1,v*b2,v*b3],
            [u*b2,v*b2,c,e],
            [u*b3,v*b3,e,g]])
F = a*v-u*u*b1
assert N*Matrix([v,-u,0,0]) == Matrix([F,0,0,0])
for i,j,R in [(2,2,b1*c-v*b2**2),(2,3,b1*e-v*b2*b3),(3,3,b1*g-v*b3**2)]:
    assert factor(N.extract([1,i],[1,j]).det()-v*R) == 0
assert N.subs({u:0,v:0,a:0})[:, :2] == Matrix.zeros(4,2)
assert [d for d in range(0,25) if d-12 >= 2 and 18-d >= 4] == [14]
print("PASS: cubic separator, rank-one residual minors, forced line and section degrees")
