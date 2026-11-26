"""Irreducible auxiliary pencils and their two nonreduced cubic inverse formulas."""
import sympy as sp
from verify_singular_integral_commutators import adjugate_map, coordinates, six_space, variables


def main():
    x, y, z = variables
    p, t, u, a, b, d = coordinates
    k = sp.symbols('k', nonzero=True)
    half = sp.Rational(1, 2)
    j = sp.Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    f = sp.Matrix([[0, 0, 0], [1, 0, 0], [0, -1, 0]])
    one = sp.Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 1]])
    def second(s):
        return sp.Matrix([[s, -s, 0], [0, -s, 1], [s*s, -s*s, 0]])

    models = [
        (j, k*f, sp.Matrix([[0, 1, 0], [0, 0, -half], [0, 0, 0]]),
         sp.Matrix([[0, 0, 0], [-k/2, 0, 0], [0, -k, 0]]),
         [-u, k*p+d, b, k*a, d, k*t], 4*k**10),
        (one, second(k), sp.Matrix([[0, 1, -1/k], [0, 1, -1/(2*k)], [0, 0, 0]]),
         sp.Matrix([[0, -k, 1], [k/2, -k, half], [0, -k**2, k]]),
         [-2*k*k*p+k*k*a-k*d+u, -k*p+k*t-b+d, p-a, -k*a+k*t-d, a-t, k*k*t],
         -4*k**10),
    ]
    for first, other, m, n, params, minor in models:
        aux = x*sp.eye(3)+y*first+z*other
        target = x*sp.eye(3)+y*m+z*n
        assert sp.factor(aux.det()-target.det()) == 0
        assert len(sp.Poly(sp.factor(sp.Matrix(params).jacobian(coordinates).det()), k).terms()) == 1
        assert sp.Matrix(params).jacobian(coordinates).subs(k, 1).det() != 0
        basis = six_space(aux)
        q = [sum((c*r[i] for c, r in zip(params, basis)), sp.zeros(3)).applyfunc(sp.factor)
             for i in range(3)]
        assert all(sp.denom(w).is_Integer for m0 in q for w in m0)
        adjugate_map(q, target, [0, 1, 2, 3, 4, 6], minor, (1, 0, 0))
        if first == j:
            change = sp.diag(1, 1, -2)
            assert change.inv()*m*change == j
            assert change.inv()*n*change == -k*f/2
        else:
            change = sp.Matrix([[1, 0, 1], [0, -1, 1], [0, -2*k, 0]])
            assert change.inv()*m*change == one
            assert (change.inv()*n*change)[0, 0] == -k/2
    print('PASS: both nonreduced adjugate identities, six-dimensional coordinates and target parameter -k/2')


if __name__ == '__main__':
    main()
