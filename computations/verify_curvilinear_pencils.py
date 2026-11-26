"""Local equations and length bounds in the marked curvilinear pencil proof."""
import sympy as sp


def check_embedded_point():
    x, y, z, w, h = sp.symbols("x y z w h")
    alpha, beta, b, c, k, kp = sp.symbols("alpha beta b c k kp")
    field = sp.QQ.frac_field(alpha, beta, b, c, k, kp)
    # Localize at the other root and at the unit solving the last row.
    unit = 1 - c*z
    simple_root = (beta-k)*unit - b*(alpha+kp*z)
    inverse = sp.Symbol("inverse")
    equations = [w-b*z-c*h, y*z-alpha*h, w*x-beta*h,
                 x-y-k*z-kp*h, h-w*z, inverse*unit*simple_root-1]
    basis = sp.groebner(equations, inverse, h, w, y, x, z, domain=field)
    for consequence in (z*z, x*z, h, w-b*z, y-x+k*z):
        assert basis.reduce(consequence)[1] == 0
    model = sp.groebner([z*z, x*z], x, z, domain=field)
    for equation in equations[:-1]:
        assert model.reduce(equation.subs({h: 0, w: b*z, y: x-k*z}))[1] == 0


def check_lengths():
    n = sp.Symbol("n", integer=True)
    hilbert = lambda d: (d+2)*(d+1)/2
    assert sp.expand(hilbert(n)-3*hilbert(n-2)+2*hilbert(n-3)) == 3
    for sites in range(1, 6):
        max_secants = sites*(sites-1)//2 + min(sites, 6-sites)
        assert max_secants <= 11
        assert (max_secants == 11) == (sites == 5)
    # D, isolated nodes, and exceptional tangent planes.
    for isolated in range(3):
        for embedded in range(3-isolated):
            capacity = 2 + 6*isolated + 3*embedded
            if capacity >= 11:
                assert (isolated, embedded) in {(1, 1), (2, 0)}
    assert (18-7, 20-7, (17+1)//2) == (11, 13, 9)


if __name__ == "__main__":
    check_embedded_point()
    check_lengths()
    print("PASS: localized embedded conic ideal, degree-three minors, and triple-point bounds")
