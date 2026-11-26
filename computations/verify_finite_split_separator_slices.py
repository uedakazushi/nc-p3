"""Check determinant and trace coefficients used in the split-separator argument."""
from sympy import Matrix, symbols

a, t, x, y, z = symbols("a t x y z")
j = Matrix(2, 2, symbols("j:4"))
b = Matrix([x, y])
star = Matrix([[a]]).row_join(b.T).col_join(b.row_join(t*j))
assert (star.det()-a*t**2*j.det()+t*(b.T*j.adjugate()*b)[0]).expand() == 0
m = Matrix(3, 3, symbols("m:9"))
n = Matrix([[0,1,0],[0,0,1],[0,0,0]])
assert ((n*n*m).trace()-m[2,0]).expand() == 0
assert ((n*n*m*m).trace()-m[2,1]*m[1,0]).subs(m[2,0], 0).expand() == 0
skew = Matrix([[0,-z,y],[z,0,-x],[-y,x,0]])
v = Matrix([x,y,z])
assert ((skew+t*m).det().expand().coeff(t,1)-(v.T*m*v)[0]).expand() == 0
print("PASS: split determinant, alternating net and nilpotent trace identities")
