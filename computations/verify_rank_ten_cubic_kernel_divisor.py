"""Exact bidegrees and restrictions of the two parametrized line families."""
import sympy as sp


def bidegree(poly, tvars, xvars):
    powers = sp.Poly(poly, *(tvars+xvars)).monoms()
    degrees = {(sum(p[:2]), sum(p[2:])) for p in powers}
    assert len(degrees) == 1
    return degrees.pop()


def main():
    s, t, x, y = sp.symbols("s t x y")
    z1 = sp.Matrix([-t*x, s*x, -t*y, s*y])
    z2 = sp.Matrix([-t*x, s*x-t*y, s*y, 0])
    assert (sp.Matrix([[s*s, s*t, t*t, 0]])*z2).applyfunc(sp.expand) == sp.zeros(1, 1)
    for z in (z1, z2):
        restriction = sum((j+1)*s**(3-j)*t**j*z[j] for j in range(4))
        assert bidegree(restriction, (s, t), (x, y)) == (4, 1)
    form = z1[0]*y - z1[2]*x
    assert bidegree(z1[1]*form-z1[3]*y*z1[0], (s, t), (x, y)) == (2, 3)
    assert bidegree(z2[0]*z2[1]*z2[2], (s, t), (x, y)) == (3, 3)
    print("PASS: limiting lines, quartic edge divisor and residual minor bidegrees")


if __name__ == "__main__":
    main()
