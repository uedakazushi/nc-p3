"""Verify the conic multiplicity and the involution orbit used for doubled binary graphs."""
from sympy import Matrix, Poly, expand, symbols

p, q, s, t, k = symbols("p q s t k")
# For q=k*p, symmetric coordinates s=p+q,t=p*q give the conic.
conic = k*s**2 - (1+k)**2*t
assert expand(conic.subs({s:p+q,t:p*q}) - (q-k*p)*(k*q-p)) == 0
pullback = Poly(conic.subs({s:p+q,t:p*q}), q)
assert pullback.diff().eval(k*p).expand() == (p*(k**2-1)).expand()
# k=-1 gives an involution and the conic specializes to a double line.
assert conic.subs(k,-1) == -s**2

a,b = symbols("a b", nonzero=True)
A,B = Matrix([[0,a],[1,0]]),Matrix([[0,b],[1,0]])
assert A*A == a*Matrix.eye(2) and B*B == b*Matrix.eye(2)
assert A*B == Matrix.diag(a,b)
# Common fixed-point case, including noncommuting involutions with parabolic product.
A,B = Matrix([[-1,1],[0,1]]),Matrix([[-1,3],[0,1]])
assert A*B == Matrix([[1,-2],[0,1]])
assert A[1,0] == B[1,0] == 0
print("PASS: doubled conic pullback, involution specialization, exchanging and fixed orbits")
