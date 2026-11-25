"""Mixed-jet and nonisomorphism checks for finite planar window invariants."""

import sympy as sp
from finite_planar_window_orbits import first_relations, flagged_key, monomials


def main():
    x, y = sp.symbols("x y")
    types = [(2, 3)] + [(0, 1)] * 16
    mons = monomials(2, 3)

    def coefficients(f):
        poly = sp.Poly(f, x, y)
        return sp.Matrix([poly.coeff_monomial(x**a * y**b) for a, b in mons])

    def operator(f, h, unit=1):
        return sp.Matrix.hstack(*(coefficients(unit * f**a * h**b) for a, b in mons))

    local = sp.Matrix.hstack(coefficients(x * x), coefficients(y * y))
    assert sp.factor(first_relations(local, mons) / (x * y)) in (4, -4)
    ideal = local.col_join(sp.zeros(16, 2))
    rep = sp.Matrix.hstack(*(coefficients(f) for f in [1, x, y, x * y]))
    coords = [
        rep.col_join(sp.Matrix([[1, t, t * t, t**3] for t in range(1 + s, 17 + s)]))
        for s in range(3)
    ]
    key = flagged_key(types, ideal, coords)
    f, h = 2 * x + 3 * y + x * y, -y + x * x
    moved_ideal = sp.diag(operator(f, h), sp.eye(16)) * ideal
    moved = []
    for s, b in enumerate(coords):
        g = sp.diag(operator(f, h, 3 + s + x - y), *range(1 + s, 17 + s))
        basis = sp.eye(4)
        basis[1, 2] = 2
        moved.append(g * b * basis + moved_ideal * sp.ones(2, 4))
    assert key == flagged_key(types, moved_ideal, moved)
    # A reparametrization at just one position changes the normalized record.
    assert key != flagged_key(
        types,
        ideal,
        [coords[0], sp.diag(operator(x, 2 * y), sp.eye(16)) * coords[1], coords[2]],
    )
    # Grade-d initial forms can be combinations of columns with cancelled lower terms.
    larger = monomials(2, 5)
    fcol = lambda f: sp.Matrix(
        [sp.Poly(f, x, y).coeff_monomial(x**a * y**b) for a, b in larger]
    )
    pencil = sp.Matrix.hstack(fcol(x**3 + y**4), fcol(y**3 + x**3))
    assert sp.factor(first_relations(pencil, larger) / (x * x * y * y)) in (9, -9)
    print(
        "PASS: planar Wronskians, common nonlinear change of parameters, unit and ideal lifts, distinct relative tangent data"
    )


if __name__ == "__main__":
    main()
