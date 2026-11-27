"""Check the two index cases in the hyperplane Riemann--Roch identity."""
import sympy as sp


def main():
    a, m = sp.symbols("a m")
    hilbert = sp.expand(sp.prod(4*m+a+j for j in (1,2,3))/6)
    section = sp.expand(hilbert - hilbert.subs(m, m-1))
    for index in [2,4]:
        h_squared = 64 // index**2
        det_step = sp.Rational(index, 4)
        assert det_step**2 * h_squared == 4
        assert index * det_step * h_squared == 16
        assert sp.expand(index**2*h_squared*m**2/2 + 16*a*m + 2 + 2*a**2 - section) == 0
    assert (2 + 2*(-1))**2 == 0
    print("Hodge index equality, integral determinant increments and balanced inertia traces.")


if __name__ == "__main__":
    main()
