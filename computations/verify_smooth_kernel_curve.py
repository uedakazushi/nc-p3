"""Check the cone genus alternative and the Koszul cubic Chern character."""
from itertools import combinations
from math import factorial
from sympy import expand, prod, symbols

b = symbols("b")
for c, target in [(2*b, -(b-2)**2), (2*b+1, -(b-1)*(b-2))]:
    assert expand(c - (b-1)*(c-b-1) - 3 - target) == 0
roots = symbols("x:6")
coefficients = (20, -10, 4, -1)
for degree in range(4):
    character = sum(coefficients[k] * (-sum(subset))**degree
                    for k in range(4) for subset in combinations(roots, k)) / factorial(degree)
    expected = sum(prod(c) for c in combinations(roots, 3)) if degree == 3 else 0
    assert expand(character - expected) == 0
print("PASS: both smooth cone classes have nonpositive defect; degree-three Koszul character equals c3")
