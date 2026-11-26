"""Polynomial identities and the Jordan-product obstruction for the rank-nine reduction."""
from sympy import Matrix, expand, symbols

a,b,c,d,e,h,k = symbols("a b c d e h k")
N = Matrix([[0,1,0],[0,0,1],[0,0,0]])
Y = Matrix([[a,b,d],[c,-a-h,e],[0,-c,h]])
assert expand((N*Y**2).trace() - c*(a-h)) == 0
Y = Y.subs(h,a)
assert expand((Y**2).trace()/2 - (3*a*a+c*(b-e))) == 0
Y = Y.subs({a:k*c,e:b+3*k*k*c})
assert expand(Y.det()-c*c*(k**3*c-d)) == 0
# Two bidirected pairs force a nonzero third pair under their Jordan product.
x,y,z,w=symbols("x y z w", nonzero=True)
A=Matrix([[0,x,0],[y,0,0],[0,0,0]])
B=Matrix([[0,0,0],[0,0,z],[0,w,0]])
assert A*B+B*A == Matrix([[0,0,x*z],[0,0,0],[y*w,0,0]])
print("PASS: nilpotent-plane polarizations, determinant, and third Jordan pair")
