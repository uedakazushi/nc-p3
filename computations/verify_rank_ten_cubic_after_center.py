"""Verify the pencil normal forms and their degree and evaluation identities."""
import sympy as sp


def coefficient_rank(columns, variables, monomials):
    return sp.Matrix([
        [sp.Poly(entry, *variables).coeff_monomial(m)
         for entry in column for m in monomials] for column in columns
    ]).rank()


def main():
    s, t, a, b, f, g = sp.symbols("s t a b f g")
    u, v0, v1 = (sp.eye(3).col(i) for i in range(3))
    v = s * v0 + t * v1
    first = sp.Matrix.hstack(s*s*u, s*t*u, s*v, t*v)
    second = sp.Matrix.hstack(s*s*u, s*t*u, t*t*u, s*v)
    for matrix in (first, second):
        cols = [matrix.col(i) for i in range(4)]
        assert matrix.rank() == 2
        assert matrix.subs({s: 0, t: 1}).rank() == 1
        assert coefficient_rank(cols, (s, t), (s*s, s*t, t*t)) == 4
        multiples = [l * col for l in (s, t) for col in cols]
        assert coefficient_rank(multiples, (s, t), (s**3, s*s*t, s*t*t, t**3)) == 6
        values = [sp.Matrix.vstack(s*col, t*col) for col in cols]
        coefficients = [value.applyfunc(lambda x: sp.Poly(x, s, t).coeff_monomial(m))
                        for value in values for m in (s**3, s*s*t, s*t*t, t**3)]
        assert sp.Matrix.hstack(*coefficients).rank() == 5
    point = sp.Matrix([-b*f, a*f, -b*g, a*g])
    assert (first * point - (t*a-s*b)*(s*f*u+g*v)).applyfunc(sp.expand) == sp.zeros(3, 1)
    print("PASS: both cubic corridor restrictions, five-dimensional image and successor projection")


if __name__ == "__main__":
    main()
