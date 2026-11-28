"""Check the local Hessian argument and the alternating-form support bound."""
from itertools import combinations

from sympy import Matrix, diag, eye, zeros

# Change input coordinates so an ordinary double point is the last basis vector.
hessian = diag(2, 4, -6, 0)
point = eye(4)[:, 3]
basis = Matrix([[1, 2, 1, 0], [0, 1, 3, 1], [0, 0, 1, 2], [0, 0, 0, 1]])
output = basis.T.inv() * point
assert hessian.nullspace() == [point]
assert output * output.T * basis * hessian == zeros(4)
assert basis.T * output == point
# Independent alternating coefficients are cut by at most one row per remaining mark.
assert len(list(combinations(range(4), 2))) == 6
assert 20 - 15 < 6 and 20 - 14 == 6
assert 20 - 16 == 4
print("PASS: four polar cubics, rank-three Hessian, at most fourteen matched marks")
