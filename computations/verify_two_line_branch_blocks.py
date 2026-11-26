"""Check the Hankel rank locus and the opposite submatrix on two lines."""
from sympy import Matrix, eye, groebner, kronecker_product, symbols

a, b, c, d = symbols("a b c d")
hankel = Matrix([[a, b, c], [b, c, d]])
minors = [hankel[:, cols].det() for cols in ((0, 1), (0, 2), (1, 2))]
s, t = symbols("s t")
param = {a: s**3, b: s**2*t, c: s*t**2, d: t**3}
assert all(f.subs(param) == 0 for f in minors)
# On the a=1 chart the ideal cuts out the graph of (b^2,b^3).
ideal = groebner([f.subs(a,1) for f in minors], c,d,b)
assert ideal.reduce(c-b**2)[1] == ideal.reduce(d-b**3)[1] == 0
assert kronecker_product(Matrix([[0,1],[-1,0]]), eye(3)).rank() == 6
print("PASS: twisted-cubic Hankel locus and nondegenerate six-row submatrix")
