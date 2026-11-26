"""The balanced inverses for the tacnode, triangle, and ordinary triple point."""
import sympy as sp
from verify_singular_integral_commutators import (
    adjugate_map, coordinates, six_space, variables,
)


def models():
    x, y, z = variables
    p, t, u, a, b, d = coordinates
    lam = sp.symbols('l')
    tangent = sp.Matrix([[y, -x, 0], [z, 0, -x], [0, y+z, -x-y]])
    tangent_target = sp.Matrix([[x+y, -x, 0], [-x+z, 2*x, -x], [0, -x+z, x-y]])
    yield tangent, tangent_target, [u, -t+d, -d, a-b, b, p], (1, 1, 0), [0, 1, 2, 3, 4, 6], -1

    triangle = sp.Matrix([[x, lam*z, 0], [0, y, x], [y, 0, z]])
    triangle_target = sp.Matrix([[x, -z, 0], [0, lam*y, -x], [y, 0, -lam*z]])
    yield triangle, triangle_target, [u, -d, (lam-1)*a, t, -b, lam*p], (1, 1, 1), [0, 1, 2, 3, 5, 6], -lam**5

    triple = sp.Matrix([[y, -x, 0], [z, x, -x], [0, x+z, -y]])
    triple_target = sp.Matrix([[y, -x, 0], [x+z, -x, -x], [0, -x+y+z, -y]])
    yield triple, triple_target, [u, -t, -d, p, b, a], (1, 2, 1), [0, 1, 2, 3, 4, 6], -1


def check_tangency():
    n, m, j, k, l, c, f = sp.symbols('n m j k l c f', nonzero=True)
    upper = -k*c/l+j*(n-m)*f/(n*l)
    lower = j*(n-m)*c/(m*k)-l*f/k
    leading = sp.factor((upper*lower-c*f)*(-m*n*k*l)/(j*(n-m)))
    assert sp.expand(leading-(n*k*c**2-j*(n-m)*c*f+m*l*f**2)) == 0
    a, b = sp.symbols('a b')
    cut = -m*k*a*a+j*(m-n)*a*b-n*l*b*b
    assert sp.expand(sp.discriminant(cut.subs(b, 1), a)-sp.discriminant(leading.subs(f, 1), c)) == 0
    x, y, z = variables
    block = sp.Matrix([[1, 2, -2], [-2, 2, 1], [-3, 2, -3]])
    off1 = sp.Matrix([[0, 1, 0], [0, 0, 0], [-1, 0, 0]])
    off2 = sp.Matrix([[0, 0, 1], [-2, 0, 0], [0, 0, 0]])
    conic = -30*x*x-15*x*y-2*x*z-(y+2*z)**2
    assert sp.expand((x*block+y*off1+z*off2).det()-x*conic) == 0
    assert sp.hessian(conic, variables).det() == 1568
    unknowns = sp.symbols('a:h')
    first, second = [sp.Matrix(2, 2, unknowns[i:i+4]) for i in (0, 4)]
    residual = sp.Matrix([[4, 1], [2, -4]])
    lin = [*(first[:, 0]*2-second[:, 1]-block[1:, 0]),
           *(second[0, :]-first[1, :]-block[0, 1:]), (first*residual).trace(),
           (second*residual).trace()]
    sol, = sp.linsolve(lin, unknowns)
    g, h = unknowns[-2:]
    parabola = (g-h)**2-10*g-4*h+11
    rem = (first*second-second*first-residual).subs(dict(zip(unknowns, sol)))
    assert (rem+residual*parabola/8).applyfunc(sp.expand) == sp.zeros(2)
    assert sp.hessian(sp.Poly(parabola, g, h).homogenize(z).as_expr(), (g, h, z)).det()


def main():
    lam = sp.symbols('l')
    x, y, z = variables
    for aux, target, params, point, rows, minor in models():
        basis = six_space(aux)
        q = [sum((c*m[i] for c, m in zip(params, basis)), sp.zeros(3)).applyfunc(sp.factor)
             for i in range(3)]
        assert all(sp.denom(h).is_Integer for m in q for h in m)
        packed = sp.Matrix.vstack(*[m.reshape(9, 1) for m in q]).jacobian(coordinates)
        assert packed.subs(lam, 1).rank() == 6
        phi = adjugate_map(q, target, rows, minor, point)
        assert all(sp.denom(h).is_Integer for h in phi)
        if lam in aux.free_symbols:
            assert sp.expand(target.subs(lam, 1).det()) == 0
    # The two constant-right, moving-left tacnodal models have a common invariant line.
    for aux in [sp.Matrix([[y, -x, 0], [z, -y, 0], [0, z-h*y, -x]]) for h in (0, 1)]+[
            sp.Matrix([[y, -x, 0], [0, y, -x], [z, 0, -x]])]:
        assert all(m[0, 2] == m[1, 2] == 0 for q in six_space(aux) for m in q)
    check_tangency()
    print('PASS: three balanced adjugate inverses, excluded tacnodal auxiliaries and affine parabola')


if __name__ == '__main__':
    main()
