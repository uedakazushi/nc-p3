"""First derivative of the elliptic multiplication pairing at a nontrivial two-torsion point."""

import numpy as np
import sympy as sp

from verify_elliptic_bezout_torsion import P, inverse_moment, multiplication, verify_quadratic_presentation


def main():
    x, y = sp.symbols("x y")
    cubic = x**3-x
    base = [1, x, y, x**2]
    sections = [1, x, y, y/x]
    # With local parameter t=y at (0,0), x(t) has no linear term.
    derivatives = [0, 0, 0, 1/x]
    symmetric = multiplication(base, base, 1, x, y, cubic)
    branch = multiplication(sections, sections, x, x, y, cubic)
    # The two nonsymmetric summands may have a pole separately. Subtract before reduction.
    octic_basis = [1, x, x**2, x**3, x**4, y, x*y, x**2*y]
    derivative = np.zeros((4, 4, 8), dtype=np.int64)
    for i in range(4):
        for j in range(4):
            value = sp.cancel(x*(sections[i]*derivatives[j]-derivatives[i]*sections[j]))
            poly = sp.Poly(value, x, y)
            derivative[i, j] = [int(poly.coeff_monomial(m)) for m in octic_basis]
    assert np.array_equal(derivative, -derivative.transpose(1, 0, 2))

    def evaluate(a, transform):
        b = transform@a
        ordinary = np.einsum("j,ijk->ik", a, branch)
        first_order = np.einsum("j,ijk->ik", b, branch)+np.einsum("j,ijk->ik", a, derivative)
        return inverse_moment(np.vstack((ordinary, first_order)) % P, symmetric)

    verify_quadratic_presentation(evaluate, 581)
    print("Verified the derivative Bezout matrix at two-torsion and its nineteen independent quartic minors.")


if __name__ == "__main__":
    main()
