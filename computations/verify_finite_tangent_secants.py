"""Check the quadratic residue bound for both splitting types."""
import sympy as sp

s, t = sp.symbols("s t")
examples = [
    [[s, 0, 0], [0, s, 0], [0, t, s], [0, 0, t]],
    [[s, 0, 0], [0, s, 0], [0, t, 0], [0, 0, s], [0, 0, t]],
]
for pencil in examples:
    rows = [[sp.sympify(x) for x in row] for row in pencil]
    columns = []
    for row in rows:
        for variable in (s, t):
            polys = [sp.Poly(variable * x, s, t) for x in row]
            columns.append([p.coeff_monomial(m) for p in polys for m in (s*s, s*t, t*t)])
    matrix = sp.Matrix(columns)
    expected = sp.zeros(9, 1)
    expected[2] = 1
    assert matrix.nullspace() == [expected]
print("PASS: both splitting types have exactly the exceptional-direction quadratic quotient")
