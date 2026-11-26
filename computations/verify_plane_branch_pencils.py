"""Check the two rulings of a Segre factor, a cone syzygy and a plane pencil."""
from sympy import Matrix, eye, symbols
x,y = symbols("x y")
def image_rank(linear):
    return Matrix.hstack(*(Matrix(row).reshape(len(row),1) for row in linear)).rank()
# Three vector-valued linear forms with a linear syzygy generate a full
# three-dimensional pure space on multiplication by a second binary letter.
forms = [Matrix([x,0,0]),Matrix([y,0,0]),Matrix([0,x,y])]
def coefficients(f):
    return [s.expand().coeff(x,2) for s in f]+[s.expand().coeff(x,1).coeff(y,1) for s in f]+[s.expand().coeff(y,2) for s in f]
assert image_rank([coefficients(a*f) for a in (x,y) for f in forms]) == 5
assert image_rank([coefficients(a*f) for a in (x,y) for f in forms[:2]]) == 3
# A smooth ruling gives an unrestricted spectator factor of dimension two.
K = Matrix([[1,2,0],[0,1,1],[1,0,2],[0,1,0]])
assert eye(2).rank()*K.rank() == 6
mat = Matrix([[x,y,0],[0,0,x],[0,0,y],[0,0,0]])
assert Matrix.vstack(mat.diff(x),mat.diff(y)).rank() == 3
assert mat.subs({x:1,y:2}).rank() == mat.subs({x:0,y:1}).rank() == 2
print("PASS: smooth spectator ranks, the pure cone subspace and the plane pencil")
