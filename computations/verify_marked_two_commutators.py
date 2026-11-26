"""Recover commuting mark coordinates for both semisimple and Jordan generators."""
from sympy import Matrix, eye

def traceless(a):
    return a-eye(3)*a.trace()/3

for x, z in (
    (Matrix.diag(1, 2, -3), Matrix([[0,1,3],[2,0,1],[1,1,0]])),
    (Matrix([[0,1,0],[0,0,1],[0,0,0]]), Matrix([[0,0,0],[1,0,0],[0,3,0]]))
):
    y = traceless(2*x*x+3*x)
    c0, c1 = x*z-z*x, y*z-z*y
    assert c0.det() != 0
    columns = []
    for i in range(9):
        unit = Matrix.zeros(3)
        unit[i] = 1
        columns.append((unit*c0+c0*unit).reshape(9, 1))
    d = (Matrix.hstack(*columns).inv()*c1.reshape(9, 1)).reshape(3, 3)
    assert traceless(d)/2 == x and traceless(d*d)/2 == y
    assert (d*d*c0).trace() == 0
print("PASS: unique Sylvester recovery for distinct and repeated eigenvalues")
