"""Check the admissible ranks for a two-line corridor and the restricted slice."""
from sympy import Matrix, eye, kronecker_product

ranks = [(a, b, r, 6+a+b-r) for a in range(2, 5) for b in range(2, 5)
         for r in range(1, min(2, a-2, b-2)+1)]
assert min(row[3] for row in ranks) == 11
assert [row[:3] for row in ranks if row[3] == 11] == [(3, 3, 1)]
bracket = Matrix([[0, 1], [-1, 0]])
assert kronecker_product(bracket, eye(3)).rank() == 6
print("PASS: full middle restrictions give a unique rank-eleven endpoint profile")
