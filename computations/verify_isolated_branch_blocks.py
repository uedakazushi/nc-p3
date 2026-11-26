"""Check maximal-minor degrees for all four rational normal successor models."""
from functools import reduce
from sympy import Matrix, gcd, symbols

x, y = symbols("x y")
for degree in range(4):
    rows = []
    for j in range(degree):
        row = [0]*4
        row[j:j+2] = [y, -x]
        rows.append(row)
    for j in range(degree + 1, 4):
        row = [0]*4
        row[j] = x + j*y
        rows.append(row)
    matrix = Matrix(rows)
    minors = [(-1)**j * matrix[:, [i for i in range(4) if i != j]].det()
              for j in range(4)]
    common = reduce(gcd, minors)
    assert common.as_poly(x, y).total_degree() == 3 - degree
    assert matrix * Matrix(minors) == Matrix.zeros(3, 1) or all(
        a.expand() == 0 for a in matrix * Matrix(minors))
    nonzero = [f for f in minors if f != 0]
    assert len(nonzero) == degree + 1
print("PASS: isolated successor pencils yield degree d and maximal-minor divisor 3-d for d=0,1,2,3")
