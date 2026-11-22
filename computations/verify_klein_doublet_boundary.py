"""Verify the completion identities, the zero pattern and the exceptional limit."""
from itertools import product
import sympy as sp
from verify_klein_two_doublets import cross_blocks

u,v,w,d,p,q,r,s = sp.symbols('u v w d p q r s')
U,V = sp.Matrix([[u,v],[w,d]]),sp.Matrix([[p,q],[r,s]])
X,Y = cross_blocks(U,V)
m,n = sp.symbols('m n')
res = Y*sp.diag(n,m)*X
coeff,_ = sp.linear_eq_to_matrix([res[0,1],res[1,0]],[m,n])
assert sp.expand(coeff.det()-X.det()*(p*s*v*w+r*q*u*d)) == 0

for bits in product(range(2),repeat=8):
    us,vs = sp.Matrix(2,2,bits[:4]),sp.Matrix(2,2,bits[4:])
    if all(bits) or any(sum(row)<1 for t in (us,vs) for z in (t,t.T) for row in z.tolist()):
        continue
    if any(sum(t[i,i^h] for t in (us,vs) for i in (0,1))<2 for h in (0,1)):
        continue
    # For zero patterns it is only the vanishing of both determinant monomials that matters.
    if not (us[0,1]*us[1,0]*vs[0,0]*vs[1,1] or us[0,0]*us[1,1]*vs[0,1]*vs[1,0]):
        assert list(us).count(0)==list(vs).count(0)==1
        i,j = next((i,j) for i,j in product(range(2),repeat=2) if us[i,j]==0)
        assert vs[j,i]==0 or vs[j^1,i^1]==0

a,b,c,d,z,h = sp.symbols('a b c d z h')
x,y = cross_blocks(sp.ones(2),sp.Matrix([[1,1],[1,z]]))
k = (z*b+c)/(z+1)
center = sp.Matrix([[k,b],[c,k]])
diag = y*center.adjugate()*x
assert sp.cancel(diag[0,1]) == sp.cancel(diag[1,0]) == 0
aa = center+sp.diag(v,w)
bb = center+h*sp.diag(v,w)
assert sp.cancel(bb.det()-h*aa.det()-(h-1)*(h*v*w-center.det())) == 0
assert sp.limit(k+(k*k-b*c)/(d-k),z,-1) == -d
assert sp.limit(k+(k*k-b*c)/(a-k),z,-1) == -a
print('Singular zero patterns, completion determinant, central quadric and both exceptional limits verified.')
