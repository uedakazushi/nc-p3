"""Matrix identities for the line and transverse conic commutator inverse."""
import sympy as sp
from verify_prescribed_discriminant_seeds import bracket
from verify_singular_integral_commutators import (
    coordinates, equal_polynomials, positions, six_space, symmetric, variables,
)


def check_affine_conic():
    x, y, z = variables
    q0 = sp.diag(2, -1, -1)/3
    q1 = sp.Matrix([[4, 1, 1], [-1, -2, 3], [0, 3, -2]])
    q2 = sp.Matrix([[-1, -6, 9], [3, 2, -3], [-6, -9, -1]])/3
    c0 = bracket(q1, q2)
    target = x*c0+y*bracket(q2, q0)+z*bracket(q0, q1)
    conic = 38*x*x-260*x*y+308*x*z+27*y*y-30*y*z-3*z*z
    assert sp.expand(target.det()-x*conic) == 0
    assert sp.hessian(conic, variables).det() and sp.discriminant(conic.subs({x: 0, z: 1}), y) == 1224
    rows = [q[0, 1:] for q in (q1, q2)]
    cols = [q[1:, 0] for q in (q1, q2)]
    assert sp.Matrix.vstack(*rows).det() and sp.Matrix.hstack(*cols).det() and c0[0, 0]
    residual = c0[1:, 1:]-cols[0]*rows[1]+cols[1]*rows[0]
    coeff = sp.symbols('a:h')
    first, second = sp.Matrix(2, 2, coeff[:4]), sp.Matrix(2, 2, coeff[4:])
    linear = [*(rows[0]*second-rows[1]*first-c0[0, 1:]),
              *(first*cols[1]-second*cols[0]-c0[1:, 0]),
              (first*residual).trace(), (second*residual).trace()]
    solution, = sp.linsolve(linear, coeff)
    g, h = coeff[-2:]
    fiber = 27*g*g-24*g*h+234*g-6*h*h-140*h+459
    remainder = (bracket(first, second)-residual).subs(dict(zip(coeff, solution)))
    assert (remainder+2*fiber*residual/(45*residual[0, 0])).applyfunc(sp.expand) == sp.zeros(2)
    homogeneous = sp.Poly(fiber, g, h).homogenize(z).as_expr()
    assert sp.hessian(homogeneous, (g, h, z)).det() and sp.discriminant(homogeneous.subs({z: 0, h: 1}), g)


def main():
    x, y, z = variables
    p, t, u, a, b, d = coordinates
    lam = sp.symbols('l')
    balanced = sp.Matrix([[-x, y, 0], [z, x*(1/lam-1), -y/lam], [0, z, -x]])
    target = sp.Matrix([[-lam*x, y, 0], [z, x*(lam**2-1)/lam, -y], [0, -lam*z, x]])
    basis = six_space(balanced)
    # These coordinates name the same intrinsic entries for every free parameter.
    sample = sp.Matrix.hstack(*[sp.Matrix([q1[1, 2], q2[0, 1], q2[0, 2],
                                           -q0[2, 1]/lam, q2[1, 2], q2[2, 1]])
                                for q0, q1, q2 in basis])
    parameters = sp.Matrix([-u, -(lam-1)*d-p, -b, a, d, (lam-1)*t])
    assert sp.factor(parameters.jacobian(coordinates).det()) in (lam-1, 1-lam)
    coef = sample.inv()*parameters
    q = [sum((c*m[i] for c, m in zip(coef, basis)), sp.zeros(3)).applyfunc(sp.factor)
         for i in range(3)]
    assert all(sp.denom(h).is_Integer for m in q for h in m)
    cs = [bracket(q[(i+1) % 3], q[(i+2) % 3]) for i in range(3)]
    multiplier = cs[0]*target.diff(x).inv()
    polys = [sp.Poly(h, *coordinates) for h in multiplier]
    cols = [sp.Matrix(3, 3, [sp.factor(k.coeff_monomial(m)) for k in polys])
            for m in (t*u, p*u, p*t, b*d, a*d, a*b)]
    adj = symmetric.adjugate()
    assert equal_polynomials(multiplier-sum((m*adj[i, j] for m, (i, j) in zip(cols, positions)),
                                           sp.zeros(3)))
    assert all(equal_polynomials(c-multiplier*target.diff(var)) for c, var in zip(cs, variables))
    flat = sp.Matrix.hstack(*[m.reshape(9, 1) for m in cols])
    assert sp.factor(flat[[0, 1, 2, 3, 4, 6], :].det()) == -lam**4
    assert all(sp.denom(v).is_Integer for v in flat)
    assert sp.factor(target.subs(lam, -1).det()) == 0
    actual = (flat*sp.Matrix(coordinates)).reshape(3, 3)
    change = dict(zip(coordinates, [p, t, u, -a, -b, 2*p-d]))
    other = (flat.subs(lam, -lam)*sp.Matrix(coordinates).xreplace(change)).reshape(3, 3)
    assert equal_polynomials(other-sp.diag(1, -1, 1)*actual*sp.diag(-1, 1, 1))
    assert equal_polynomials(sp.Matrix([4*actual.det()-lam*(lam+1)**2*symmetric.det()
                                        +lam*(lam-1)**2*symmetric.xreplace(change).det()]))
    # Constant right kernel on the line, with moving left kernel.
    for aux, i in [(sp.Matrix([[-x, y, 0], [z, -x, 0], [lam*y, z, -x]]), 2),
                   (sp.Matrix([[-x, y, 0], [0, -x, y], [0, z, -x]]), 0)]:
        basis = six_space(aux)
        columns = []
        for triple in basis:
            for m in triple:
                assert all(m[j, i] == 0 for j in range(3) if j != i)
                for entry in m:
                    assert len(sp.Poly(sp.denom(sp.cancel(entry)), lam).terms()) == 1
            columns.append(sp.Matrix.vstack(*[m.reshape(9, 1) for m in triple]))
        flat = sp.Matrix.hstack(*columns)
        # The nullspace free coordinates have a constant unit minor (use sl3 coordinates).
        extract = [i for i in range(flat.rows) if sum(v != 0 for v in flat[i, :]) == 1
                   and any(v in (1, -1, sp.Rational(2, 3), sp.Rational(-1, 3))
                           for v in flat[i, :])]
        assert flat[extract, :].rank() == 6
    check_affine_conic()
    print('PASS: balanced conic sheets, determinant pencil and both excluded auxiliary types')


if __name__ == '__main__':
    main()
