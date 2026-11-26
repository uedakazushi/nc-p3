"""Local rank-one ideal and the restriction of the common rows to two planes."""
from itertools import combinations

import sympy as sp

a, x, y, t, l, u, v = sp.symbols("a x y t l u v")
pencil = sp.Matrix([[a, x, y], [x, t, t], [y, 0, t]])
pairs = list(combinations(range(3), 2))
minors = [pencil.extract(i, j).det() for i in pairs for j in pairs]
ideal = sp.groebner([m.subs(a, 1) for m in minors], x, y, t)
assert ideal.polys == sp.groebner([x*x, x*y, y*y, t], x, y, t).polys
for var in (x, y, t):
    others = [z for z in (x, y, t) if z != var]
    assert sp.groebner([m.subs({a: 0, var: 1}) for m in minors], *others).contains(sp.Integer(1))
common = (sp.Matrix([[l, u, v]]) * pencil).subs({l: 0, t: 0})
assert common == sp.Matrix([[u*x + v*y, 0, 0]])
print("PASS: square-zero triple scheme and the one surviving shared plane equation")
