"""Symbolic identities used in the four-dimensional marked-pair proof."""
from itertools import combinations

import sympy as sp


def bracket(x, y):
    return x * y - y * x


def unit(row, column):
    matrix = sp.zeros(3)
    matrix[row, column] = 1
    return matrix


j = sp.Matrix(2, 2, sp.symbols("j0:4"))
basis = [sp.diag(1, 0), sp.Matrix([[0, 1], [1, 0]]), sp.diag(0, 1)]
u, v, w = [entry * j for entry in basis]
value = sp.factor((u * bracket(v, w)).trace())
assert sp.simplify(value - (j[1, 0] - j[0, 1]) * j.det()) == 0 or sp.simplify(
    value + (j[1, 0] - j[0, 1]) * j.det()
) == 0
flag_basis = [unit(0, 1), unit(0, 2) + unit(1, 1), unit(1, 2)]
assert [bracket(x, y) for x, y in combinations(flag_basis, 2)] == [
    unit(0, 1), unit(0, 2), unit(1, 2)
]
p, q, s, t = sp.symbols("p q s t")
a = p * unit(0, 1) + q * unit(1, 0)
b = s * unit(1, 2) + t * unit(2, 1)
assert bracket(a, b) == p * s * unit(0, 2) - q * t * unit(2, 0)
print("PASS: Cartan form, flag brackets, and oriented triangle identity")
