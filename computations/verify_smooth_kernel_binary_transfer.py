"""Verify residue identities, invariant-factor bounds, edge brackets and excess lengths."""
from sympy import Matrix, Poly, expand, rem, symbols, zeros

t,x,y,z,k = symbols("t x y z k")
def top(f, p):
    return Poly(rem(f,p,t),t).nth(Poly(p,t).degree()-1)
for p,q in [((t-1)**5,t-1),((t+2)**4,(t+2)**2),((t-2)**3,(t-2)**3)]:
    m,e = Poly(p,t).degree(),Poly(q,t).degree()
    for j in range(12):
        assert top(Poly(p,t).exquo(Poly(q,t)).as_expr()*t**j,p) == top(t**j,q)
    if m >= 4:
        assert 3*(m-4) < m-e
p = (t-2)**4*(t+1)**2
pairing = Matrix(3,6,lambda i,j:top(t**(i+j),p))
assert pairing.rank() == 3 and pairing[:,:3] == zeros(3)
degree = Poly((x+y+z)**5-k*x*(x+y+z)**4,x,y,z)
assert degree.coeff_monomial(x*y*z**3) == 20-4*k

# Exponential cross case: one B bracket at each F edge, one A bracket at
# each G edge. No two brackets consume the same binary factor of an even slot.
for length in range(3,26):
    degrees = [[0,0] if i%2 == 0 else [0] for i in range(length)]
    for i in range(length-1):
        even,odd = (i,i+1) if i%2 == 0 else (i+1,i)
        which = 1 if odd%4 == 1 else 0
        degrees[even][which] += 1
        degrees[odd][0] += 1
    assert all(degrees[i] == ([1,1] if i%2 == 0 else [2]) for i in range(1,length-1))
print("PASS: residue duals at multiple roots, minor degrees, excess Chern length, exponential edge degrees")
