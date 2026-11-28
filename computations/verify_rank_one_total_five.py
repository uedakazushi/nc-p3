"""Check the numerical capacities and the four separators on quadrics through a line."""
from itertools import combinations_with_replacement
from sympy import Matrix

assert 2*15 + 2*7 + 15 == 59 < 3*20
assert 6-2 > 3
monomials = [(i, j) for i, j in combinations_with_replacement(range(4), 2) if j >= 2]
assert len(monomials) == 7
# Sample four points on two different planes through x2=x3=0.
points = ((1, 2, 1, 0), (3, 2, 1, 0), (2, 0, 0, 1), (4, 3, 0, 1))
evaluation = Matrix([[p[i]*p[j] for i, j in monomials] for p in points])
assert evaluation.rank() == 4 and 7-4 < 4
print("PASS: repeated quadric budget and four independent conditions off the base line")
