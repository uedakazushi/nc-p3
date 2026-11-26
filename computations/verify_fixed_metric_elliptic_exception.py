"""A smooth cubic and an invertible metric on its second Veronese surface."""
import sympy as sp


def main():
    x, y, z = variables = sp.symbols('x y z')
    v = sp.Matrix(variables)

    def pencil(x, y, z):
        return sp.Matrix([[x, 2*z, 2*y], [2*z, y, 2*x], [2*y, 2*x, z]])

    m = pencil(*variables)
    f = 17*x*y*z-4*(x**3+y**3+z**3)
    assert sp.expand(m.det()-f) == 0 and 17**3 != 27*4**3
    other = sp.Matrix(sp.symbols('u v w'))
    assert (pencil(*other)*v-m*other).applyfunc(sp.expand).is_zero_matrix
    adj = m.adjugate()
    polys = [sp.Poly(adj[i, j], *variables) for i, j in [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]]
    mons = [x*x, y*y, z*z, y*z, x*z, x*y]
    assert sp.Matrix([[p.coeff_monomial(q) for q in mons] for p in polys]).det() == -2744
    assert adj.subs({x: 1, y: 0, z: 0}).det() == 16
    print("PASS: smooth cubic, reciprocal kernel pencil, Veronese determinant -2744, metric determinant 16")


if __name__ == '__main__':
    main()
