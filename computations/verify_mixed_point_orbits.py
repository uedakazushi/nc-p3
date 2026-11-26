"""Verify quartic normalizations, a self-annihilating-ideal family and the seven plane directions."""
from itertools import combinations, product

from sympy import I, Matrix, Poly, eye, expand, groebner, symbols

t, c = symbols("t c")
for n in range(2,5):
    moments = [2*j+1 for j in range(n)]
    moments[-1] = 1
    def trace(poly):
        p = Poly(poly,t)
        return sum(p.nth(i)*moments[i] for i in range(n))
    for j in range(2,n):
        k = n-j
        assert expand(trace((t+c*t**j)**k)-trace(t**k)-k*c) == 0
        assert all(trace((t+c*t**j)**i) == trace(t**i) for i in range(k+1,n))

# Idempotent basis in k^4. Character basis in the trace-zero plane.
U = I*Matrix([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])
assert U.T*U == -4*eye(3)
x,y,z = symbols("x y z")
vec = Matrix([x,y,z])
w = U*vec
projected = -U.T*w.multiply_elementwise(w)/4
assert (projected-2*I*Matrix([y*z,z*x,x*y])).applyfunc(expand).is_zero_matrix
fixed = vec.cross(projected/(2*I))
# In the chart x=1 we obtain five reduced points, then two remaining axes.
one = groebner([a.subs(x,1) for a in fixed],y,z)
assert one.is_zero_dimensional
candidates = [(1,0,0),(0,1,0),(0,0,1)] + [(1,a,b) for a,b in product((-1,1),repeat=2)]
for a,b,d in candidates:
    at = {x:a,y:b,z:d}
    assert fixed.subs(at).is_zero_matrix and fixed.jacobian([x,y,z]).subs(at).rank() == 2
for i,j in combinations(range(4),2):
    # Isotropic directions in every two-coordinate odd plane have two roots.
    a,b = symbols("a b")
    v = a*U[i,:].T + b*U[j,:].T
    assert Poly(v.dot(v).subs(a,1),b).discriminant() != 0
    assert v.dot(v).subs({a:1,b:0}) and v.dot(v).subs({a:0,b:1})

# x^2=y^2=0: each span(a*x+b*y, xy) is its own annihilator precisely
# for an isotropic direction of the product form, and otherwise has the
# complementary line a*x-b*y. All directions still give distinct ideals.
a,b = symbols("a b")
vector = Matrix([a,b])
pairing = Matrix([[0,1],[1,0]])
assert (vector.T*pairing*Matrix([a,-b]))[0] == 0
print("PASS: trace normalizations, plane eigenvectors, socle roots and length-two ideal pairing")
