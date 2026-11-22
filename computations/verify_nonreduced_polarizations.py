"""Verify torus characters and an Artin reduction for the nonreduced quartics."""
from sympy import Poly, QQ, resultant, symbols

h, x, y, z = variables = symbols("h x y z")
# Four coordinates in the manuscript may be permuted. The last entry is an
# unused coordinate when the two equations involve only three variables.
models = [
    (h**2, x*y + z**2, (1, 0, 0, 0)),
    (h**2, x*y, (0, 0, 0, 1)),
    (h**2, x*y + h*z, (0, 1, 1, 2)),
    (h**2, x**2, (0, 0, 1, 0)),
    (h**2, x**2 + h*y, (0, 0, 0, 1)),
    (h*x, y**2 + h*z, (0, 0, 1, 2)),
    (h*x, y**2 + h**2, (0, 0, 0, 1)),
    (h*x, y**2, (0, 0, 0, 1)),
    (h*x, h*y + x*z, (0, 0, 1, 1)),
    (h*x, y*(h+x), (0, 0, 0, 1)),
    (h*x, h*y + x**2, (0, 0, 0, 1)),
]


def weight(polynomial, weights):
    monomials = Poly(polynomial, *variables).monoms()
    values = {sum(a*b for a, b in zip(m, weights)) for m in monomials}
    assert len(values) == 1
    return values.pop()


def main():
    u, v = symbols("u v")
    line = dict(zip(variables, (u, v, u+2*v, 3*u+v)))
    for first, second, weights in models:
        character = sum(weights) - weight(first, weights) - weight(second, weights)
        assert character != 0
        # No projective common zero on the test line, including [1:0].
        f, g = (Poly(q.subs(line), u, v, domain=QQ) for q in (first, second))
        assert f.eval({u: 1, v: 0}) or g.eval({u: 1, v: 0})
        assert resultant(f.as_expr().subs(v, 1), g.as_expr().subs(v, 1), u)
        print(first, second, "character", character)
    print("PASS: homogeneity, nonzero residue characters and proper line cuts")


if __name__ == "__main__":
    main()
