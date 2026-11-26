"""Exact cusp, irreducibility and cone Hilbert polynomial checks."""
from sympy import Matrix, Poly, diff, discriminant, expand, factor, gcd, symbols


def main():
    n, x, z, h = symbols("n x z h")
    assert expand((n+1)**2-(n-2)*(n-3)) == 7*n-5
    equation = (1-x*x)*z*z+(2-2*x**4-8*x*x)*z+(1-x*x)**3
    assert expand(equation-((1+z)**2-x*x*(x*x+z)**2-x*x*(3-3*x*x+8*z))) == 0
    disc = 16*x*x*(1+x*x)*(3*x*x-1)
    assert expand(discriminant(equation, z)-disc) == 0
    assert gcd(x*x+1, 3*x*x-1) == 1
    assert diff(equation, z).subs({x: 1, z: 0}) == -8
    assert factor(equation.subs(z, 0)) == -(x-1)**3*(x+1)**3
    assert Poly(equation.subs(z, h*h), h).primitive()[0] == 1
    r, s, u, v = symbols("r s u v")
    eq = Matrix([r*u-s*s-h*h, s*v-r*r-h*h, v*v-u*u-3*r*r+3*s*s-8*h*h])
    for sign in (-1, 1):
        at = {r: 0, s: 0, h: 0, u: 1, v: sign}
        assert eq.subs(at) == Matrix.zeros(3, 1)
        assert eq.jacobian([r, s, h, v]).subs(at).rank() == 3
    print("PASS: septic cone class, exact nonsquare discriminant, cusp order, both projection centers")


if __name__ == "__main__":
    main()
