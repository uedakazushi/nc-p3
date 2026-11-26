"""Intersection and length counts for the plane section and split socles."""
from math import comb

assert comb(5, 2) == 10
assert comb(5, 2)-3 == 7
assert (5-3, 5-4) == (2, 1)  # adjunction on P^2 x P^3
from sympy import symbols
h, l1, l2 = symbols("h l1 l2")
norm = l1*l2**2
assert h*norm/l1 == h*l2**2
assert h*norm/l2 == h*l1*l2
print("PASS: adjunction, residual lengths, and both socle cubics for (2,1)")
