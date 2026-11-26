"""Exact three-dimensional linear identities for a marked finite rank-eight window."""
from sympy import Matrix, diag, eye, expand, symbols, zeros

def cross(v):
    a,b,c = v
    return Matrix([[0,-c,b],[c,0,-a],[-b,a,0]])

def sym(i,j):
    a = zeros(3)
    a[i,j] = a[j,i] = 1
    return a

basis = [diag(1,0,-1),diag(0,1,-1),sym(0,1),sym(0,2),sym(1,2)]
unit = [eye(3).col(i) for i in range(3)]
K = [cross(v) for v in unit]
pairs = [(1,2),(2,0),(0,1)]
cols = []
for index in range(3):
    for mat in basis:
        S = [mat if i == index else zeros(3) for i in range(3)]
        entries = []
        for i,j in pairs:
            A = K[i]*S[j]-S[j]*K[i] - (K[j]*S[i]-S[i]*K[j])
            entries.extend([A[0,0],A[1,1],A[0,1],A[0,2],A[1,2]])
        cols.append(Matrix(entries))
obstruction = Matrix.hstack(*cols)
assert obstruction.det() == -3456

# Symmetrizer determinant on a diagonal first matrix. Conjugation and density
# extend its equality up to sign with det([X,Z]) to every pair of matrices.
x = diag(*symbols("x0:3"))
z = Matrix(3,3,symbols("z0:9"))
all_sym = [sym(i,i) for i in range(3)] + [sym(*p) for p in pairs]
def skew_coordinates(a):
    return [a[i,j] for i,j in pairs]
mat = Matrix.hstack(*[Matrix(skew_coordinates(h*x-x.T*h)+skew_coordinates(h*z-z.T*h)) for h in all_sym])
determinant = mat.det()
commutator = (x*z-z*x).det()
assert expand(determinant - commutator) == 0 or expand(determinant + commutator) == 0

# Skew triples have the opposite sign in the trace identity.
l1,l2,l3 = symbols("l1:4", nonzero=True)
L = diag(l1,l2,l3)
G = [cross(L*v) for v in unit]
D = Matrix(3,3,lambda i,j:(G[pairs[j][0]]*G[pairs[j][1]]-G[pairs[j][1]]*G[pairs[j][0]])[pairs[i]])
for h in basis:
    P = Matrix(3,3,lambda i,j:(h*G[i]*G[j]).trace())
    assert (D*P*D.T + D.det()*h).applyfunc(expand) == zeros(3)
print("PASS: spin-one commutators, simultaneous symmetrizer, and skew trace sign")
