"""Exact inverse and interpolation in the length-four Artin cases."""
from sympy import Matrix, Poly, eye, symbols, zeros

s, x, y, z, t = symbols("s x y z t")
for size in (1, 2, 3, 4):
    shift = zeros(size)
    for i in range(size-1):
        shift[i+1, i] = 1
    n = x*shift + y*shift**2 + z*shift**3
    numerator = sum(((-n)**i * s**(size-1-i) for i in range(size)), zeros(size))
    assert ((s*eye(size)+n) * numerator - s**size*eye(size)).applyfunc(
        lambda p: p.expand()).is_zero_matrix
# Basis 1,X,Y,XY in the noncurvilinear Gorenstein algebra X^2=Y^2=0.
n = Matrix([[0, 0, 0, 0], [x, 0, 0, 0], [y, 0, 0, 0], [z, y, x, 0]])
assert (n**2)[3, 0] == 2*x*y and (n**3).is_zero_matrix
assert ((s*eye(4)+n)*(s**2*eye(4)-s*n+n**2)-s**3*eye(4)).applyfunc(
    lambda p: p.expand()).is_zero_matrix
for poly in (t**4, t**3*(t-1), t**2*(t-1)**2, t**2*(t-1)*(t-2),
             t*(t-1)*(t-2)*(t-3)):
    p = Poly(poly, t)
    companion = zeros(4)
    for j in range(3):
        companion[j, j+1] = 1
    companion[3, :] = Matrix([[-p.nth(j) for j in range(4)]])
    e = Matrix([1, t, t**2, t**3])
    assert (t*e-companion*e).applyfunc(lambda f: Poly(f, t).rem(p).as_expr()).is_zero_matrix
    for root in p.ground_roots():
        assert (companion-root*eye(4)).rank() == 3
print("PASS: Artin inverse denominators and all length-four divisor partitions")
