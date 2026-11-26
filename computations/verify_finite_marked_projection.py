"""Check the cubic local bound and a smooth determinantal projection system."""
import sympy as sp

u, v = sp.symbols("u v")
matrix = sp.Matrix([[u, 0, 0], [v, u, 0], [0, v, u], [0, 0, v]])


def row_span(d):
    columns = []
    for i in range(d):
        for row in matrix.tolist():
            forms = [sp.Poly(u**i * v**(d - 1 - i) * f, u, v) for f in row]
            columns.append([
                f.coeff_monomial(u**j * v**(d - j))
                for f in forms for j in range(d + 1)
            ])
    return sp.Matrix(columns).rank()


assert row_span(2) == 8 and row_span(3) == 12

x, y, z, t = variables = sp.symbols("x y z t")
pencil = sp.Matrix([
    [x + t, y - z, z + 2*t],
    [y + t, z - t, x + y],
    [z + x, x - t, y + 2*z],
])
cubic = pencil.det()
grad = [cubic.diff(var) for var in variables]
for var in variables:
    affine_variables = [a for a in variables if a != var]
    ideal = sp.groebner([f.subs(var, 1) for f in grad], *affine_variables)
    assert ideal.contains(sp.Integer(1))
print("PASS: sharp cubic row bound and smooth determinantal cubic on all four charts")
