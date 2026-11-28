"""Verify projection root orders, intersection budgets and polar row elimination."""
from fractions import Fraction
from math import prod
from sympy import Matrix

for k in (1, 2, 3):
    assert all(k*(3-j) > 2 for j in range(k))
    for r in range(1, 5):
        local = r + Fraction(max(k-r, 0)**2, 2)
        assert k <= local + (Fraction(1, 2) if k-r == 1 else 0)
assert [3*b+(4-b)*(5-b)//2 for b in range(4)] == [10, 9, 9, 10]
assert 6*4//3 == 8
assert Matrix([[1,1,1], [3,2,1], [6,2,0]]).det() != 0
assert prod((3,4,5))+12*prod((3,2,1)) > prod((6,5,4))
print("PASS: triple-direction orders, common-factor budgets and local polar multiplicities")
