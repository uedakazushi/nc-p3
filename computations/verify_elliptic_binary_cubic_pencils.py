"""Check the additive identities for the secant sets, including special automorphisms."""
from itertools import product

from sympy import Matrix, symbols

# An alternating 3-by-3 pencil cannot contain an invertible matrix.
x, y, z = symbols("x y z")
assert Matrix([[0,x,y],[-x,0,z],[-y,-z,0]]).det() == 0

# Use torsion coordinates in a two-dimensional group. The sign quotient is
# the double cover p, and a point in its symmetric square is an unordered pair.
N = 24
def plus(a, b):
    return tuple((u + v) % N for u, v in zip(a, b))
def minus(a):
    return tuple(-u % N for u in a)
def times(n, a):
    return tuple(n*u % N for u in a)
def sign_class(a):
    return min(a, minus(a))
def pair(y, d):
    return tuple(sorted((sign_class(plus(y,d)), sign_class(plus(y,minus(d))))))
def half(a):
    assert all(v % 2 == 0 for v in a)
    return tuple(v//2 for v in a)

two = list(product((0, N//2), repeat=2))
def check(delta, step, alpha):
    c = [(0,0), delta]
    for j in range(2, 14):
        c.append(plus(c[j-2], step if j % 2 == 0 else minus(step)))
    d = {j:half(plus(c[j],minus(c[j-1]))) for j in range(1,len(c))}
    m = {j:plus(c[j-1],d[j]) for j in d}
    for j in range(3,10):
        assert times(2,plus(m[j],minus(m[j-2]))) == (0,0)
        assert times(2,d[j]) == times(2,d[j+4])
        assert sign_class(times(2,alpha(d[j]))) == sign_class(times(2,d[j]))
    for j in range(3,10):
        target = {pair(plus(d[j+2], e), d[j+1]) for e in two}
        for b in two:
            offset = plus(b,alpha(plus(m[j-2],minus(m[j]))))
            images = {pair(plus(alpha(plus(d[j+1],e)),offset), alpha(d[j-2])) for e in two}
            assert images == target
check((2,4), (6,0), lambda v:v)
check((4,2), (0,6), minus)
check((12,12), (0,0), lambda v:(-v[1] % N, v[0]))
print("PASS: alternating plane, two-periodic sums, four-periodic gaps, and secant-set propagation")
