"""Check the loop equations for independent, proportional, and zero block derivatives."""
from itertools import product

from sympy import Matrix, zeros

covectors = [(0, 0), (1, 0), (0, 1), (1, 1), (2, -1)]
for left, right in product(covectors, repeat=2):
    matrix = Matrix([(*right, 0, 0), (*left, *right), (0, 0, *left)])
    solution = matrix.nullspace()[0]
    assert solution != zeros(4, 1) and matrix * solution == zeros(3, 1)
    assert any(solution[:2]) or any(solution[2:])
print("PASS: a nonconstant dual-number loop exists in all tested derivative ranks")
