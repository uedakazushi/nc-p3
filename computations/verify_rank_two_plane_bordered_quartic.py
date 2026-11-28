"""Schur formulas and length budgets for the bordered quartic."""
from math import comb
from sympy import Matrix, diag, factor, symbols

h, t = symbols("h t")
C = Matrix([[2, 1, 3], [1, 0, 1], [3, 1, 1]])
w = Matrix([1, t, 1-t])
border = C.row_join(w).col_join(w.T.row_join(Matrix([[-h]])))
assert factor(h**2*border.det()+(w*w.T+h*C).det()) == 0
B = diag(1, 1, 0)
m = Matrix(symbols("a b c"))
assert (t*B.det()-(m.T*B.adjugate()*m)[0]).expand() == -m[2]**2
assert comb(6,3)//2 == 10
for d in range(4):
    assert 10-d+2*d < 14
    assert 10-d+2*d+1 < 17
    assert 6-d+3*d < 14
assert 2*6+1 < 14 and 4*6//2 < 14
assert 20-2*4+1 < 14 and 2*4+2+1 < 14
print("PASS: bordered determinant, rank-two directional square, and all reduced quartic budgets")

x, y, z = symbols("x y z")
quadric = x*y+z**2+h**2
quartic = h**2*quadric
assert factor(quartic.diff(x)) == h**2*y
assert factor(quartic.diff(y)) == h**2*x
