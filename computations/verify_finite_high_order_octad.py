"""Exact binary identities for the local edge and its linked octad."""
from sympy import Matrix, eye, symbols

x = Matrix(symbols("x:2"))
d = Matrix(symbols("d:2"))
c, k, r = symbols("c k r")
row = Matrix([[-x[1], x[0]]])
y = c*x+d
g = (row*y)[0].expand()
q = c*r+k*g
m = k*x*row-r*eye(2)
assert (m*d - (q*x-r*y)).applyfunc(lambda f: f.expand()).is_zero_matrix
assert (m.det()-r**2).expand() == 0
assert (q*q-c**2*r**2-k*g*(q+c*r)).expand() == 0
print("PASS: cubic reduction, linkage determinant and common normal square")
