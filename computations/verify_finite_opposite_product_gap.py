"""Check the cofactor and alternating-pair identities in the product-gap proof."""
import sympy as sp


def cross(v):
    a, b, c = v
    return sp.Matrix([[0, -c, b], [c, 0, -a], [-b, a, 0]])


v = sp.Matrix(sp.symbols("x y z"))
a, b, c, d, e, f = sp.symbols("a b c d e f")
metric = sp.Matrix([[a, b, c], [b, d, e], [c, e, f]])
exterior = cross(v)
assert exterior.adjugate() == v * v.T
# adj(ABC)=adj(C)adj(B)adj(A) works also for singular factors.
cofactors = exterior.adjugate() * metric.adjugate() * exterior.T.adjugate()
assert (cofactors - (v.T * metric.adjugate() * v)[0] * v * v.T).applyfunc(sp.expand).is_zero_matrix

j = sp.Matrix([[0, 1], [-1, 0]])
one = sp.diag(j, j)
two = sp.diag(2*j, 3*j)
extra = sp.diag(j, sp.zeros(2))
assert extra * (two-2*one) == sp.zeros(4)
syzygy = extra * one
x = sp.Matrix(sp.symbols("x0:4"))
pencil = ((one*x) * (two*x).T + (two*x) * (one*x).T)/2
assert syzygy.rank() == 2
assert (pencil * syzygy * x).applyfunc(sp.expand).is_zero_matrix
print("PASS: symbolic ternary adjugate and rank-two mixed syzygy")
