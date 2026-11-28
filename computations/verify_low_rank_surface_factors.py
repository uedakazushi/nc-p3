"""Check the numerical steps and the symmetric pencil for a rank-two plane."""
from itertools import product
from math import comb

from sympy import Matrix, diag, factor, symbols

splits = []
for i, j in product(range(3), repeat=2):
    dimensions = ((1+i)*(1+j), (3-i)*(3-j))
    if min(dimensions) >= 3 and max(dimensions) >= 4:
        splits.append((i, j))
assert splits == [(1, 1)]
assert 20 - 14 > 4 and 20 - 4 > 2*7
chi_sextic = 1 + comb(5, 3)
base_genus = 1 - (chi_sextic - 18)
plane_genus = comb(5, 2)
assert base_genus == 8 and (2*plane_genus-2)//(2*base_genus-2) == 1

# Verify the sign and scaling of the linear determinantal presentation.
h = symbols("h")
a, b = Matrix([1, 3, 2, 5]), Matrix([2, 1, -1, 3])
c = diag(1, 2, 3, 4)
block = c.row_join(a).row_join(b).col_join(
    a.T.row_join(Matrix([[0, -h]]))).col_join(b.T.row_join(Matrix([[-h, 0]])))
assert factor((h*c+a*b.T+b*a.T).det() + h**2 * block.det()) == 0
print("PASS: quadric factor degrees, marked bounds, ruling degree and block determinant")
