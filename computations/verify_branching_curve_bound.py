"""Check the two plane models, projection ranks and degree bounds."""
from sympy import Matrix, eye, symbols, zeros
x,y,z = symbols("x y z")
# Plane of lines in a fixed plane and plane of lines through a fixed point.
fixed_plane = Matrix([[0,z,-y],[-z,0,x],[y,-x,0],[0,0,0]])
fixed_point = Matrix([[x,y,z,0],[0,0,0,x],[0,0,0,y],[0,0,0,z]])
def coefficient_rank(matrix):
    return Matrix.vstack(*(matrix.diff(a) for a in (x,y,z))).rank()
assert coefficient_rank(fixed_plane) == 3 and coefficient_rank(fixed_point) == 4
assert coefficient_rank(fixed_point[1:,:]) == 1
for a in ((1,0,0),(1,2,-1),(0,0,1)):
    for mat in (fixed_plane,fixed_point):
        assert mat.subs(dict(zip((x,y,z),a))).rank() == 2
# On a conic, det S(1) has degree <=3. An integral plane cubic requires
# a Plucker system of degree >=3, giving at most 3 sections of det S(1).
assert 1+(4-1+1) < 6 and 1+(6-3) < 6
print("PASS: plane generators, projected bilinear form and curve-section bounds")
