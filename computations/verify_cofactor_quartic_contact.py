"""Plane-section numerics and the reconstruction of a linear map from two planes."""
from sympy import Matrix
degree_choices = [d for d in range(33) if d-12 >= 3 and 2*d <= 32]
assert degree_choices == [15,16]
assert all((24-d)//2 == 4 for d in degree_choices)
A = Matrix([[1,3,0,-2],[2,0,1,4],[-1,1,3,0]])
planes = [(0,1,2),(0,1,3)]
restrictions = [A[:,plane] for plane in planes]
second_scaled = 7*restrictions[1]
common_ratio = second_scaled[0,0]/restrictions[0][0,0]
assert second_scaled[:, :2]/common_ratio == restrictions[0][:, :2]
assert Matrix.hstack(restrictions[0],second_scaled[:,2]/common_ratio) == A
assert 4*8 == 2*16 and 4*4-15 == 1
print("PASS: sectional degrees, residual line and linear gluing on intersecting planes")
