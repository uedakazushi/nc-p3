"""Ordered commutator fibers over the node and cusp, in characteristic zero."""
from itertools import combinations_with_replacement

import sympy as sp
from verify_prescribed_discriminant_seeds import axial, bracket

x, y, z = variables = sp.symbols('x y z')
p, t, u, a, b, d = coordinates = sp.symbols('P T U A B D')
symmetric = sp.Matrix([[p, a, b], [a, t, d], [b, d, u]])
positions = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


def node(lam):
    return sp.Matrix([[0, -z, y], [-y, x, (lam-1)*z], [-z, -(lam-1)*y, lam*x]])


def cusp(free):
    if free:
        return sp.Matrix([[y, -x, -x], [z, -y, -x], [0, z, -y]])
    return sp.Matrix([[y, 0, -x], [z, -y, 0], [0, z, -y]])


def six_space(matrix):
    identity = sp.eye(3)
    elementary = [identity[:, i]*identity[:, j].T-identity*sp.Rational(i == j, 3)
                  for i in range(3) for j in range(3) if (i, j) != (2, 2)]
    equations = sp.Matrix.hstack(*[matrix.diff(var) for var in variables])
    columns = []
    for index in range(3):
        for m in elementary:
            slices = []
            for i, j in combinations_with_replacement(range(3), 2):
                h = sp.zeros(3)
                h[i, j] = h[j, i] = 1
                slices.append(equations[:, 3*index:3*index+3]*axial(h*m-m.T*h))
            columns.append(sp.Matrix.vstack(*slices))
    kernel = sp.Matrix.hstack(*columns).nullspace()
    assert len(kernel) == 6
    return [[sum((row[8*i+j]*m for j, m in enumerate(elementary)), sp.zeros(3))
             for i in range(3)] for row in kernel]


def normalized(matrix, lam, parameters, nodal):
    basis = six_space(matrix)
    columns = []
    for _, q1, q2 in basis:
        columns.append(sp.Matrix([q1[0, 2], q2[0, 1], q2[0, 2], -lam*q1[2, 1],
                                  q2[1, 2], q1[1, 2] if nodal else q2[1, 0]]))
    transform = sp.Matrix.hstack(*columns).inv()*sp.Matrix(parameters)
    triple = [sum((w*m[i] for w, m in zip(transform, basis)), sp.zeros(3)).applyfunc(sp.factor)
              for i in range(3)]
    assert all(sp.denom(entry).is_Integer for m in triple for entry in m)
    assert sp.factor(sp.Matrix(parameters).jacobian(coordinates).det()) != 0
    return triple


def equal_polynomials(matrix):
    return all(all(sp.cancel(c) == 0 for c in sp.Poly(w, *coordinates).coeffs())
               for w in matrix)


def adjugate_map(q, target, rows, expected_minor):
    commutators = [bracket(q[(i+1) % 3], q[(i+2) % 3]) for i in range(3)]
    multiplier = commutators[1]*target.diff(y).inv()
    # Each monomial below occurs in precisely one adjugate entry, with coefficient one.
    monomials = [t*u, p*u, p*t, b*d, a*d, a*b]
    polynomials = [sp.Poly(w, *coordinates) for w in multiplier]
    columns = [sp.Matrix(3, 3, [sp.factor(w.coeff_monomial(h)) for w in polynomials])
               for h in monomials]
    adj = symmetric.adjugate()
    assert equal_polynomials(multiplier-sum((m*adj[i, j] for m, (i, j) in zip(columns, positions)),
                                           sp.zeros(3)))
    for c, var in zip(commutators, variables):
        assert equal_polynomials(c-multiplier*target.diff(var))
    result = sp.Matrix.hstack(*[m.reshape(9, 1) for m in columns])
    assert sp.factor(result[rows, :].det()-expected_minor) == 0
    return result


def main():
    lam = sp.symbols('l')
    params = [-lam*((lam-1)*a+u), -(lam-1)*b+lam*t, (lam-1)**2*p+lam*d,
              -lam*a, -lam*b, (lam-1)*p]
    assert abs(sp.Poly(sp.Matrix(params).jacobian(coordinates).det(), lam).LC()) == 1
    assert sp.factor(sp.Matrix(params).jacobian(coordinates).det()/(lam**5*(lam-1))) in (1, -1)
    q = normalized(node(lam), lam, params, True)
    target = node(1/lam**2)*sp.diag(1, lam, lam**2)
    assert target.subs(lam, -1).det().expand() == 0
    phi = adjugate_map(q, target, list(range(6)), lam**2*(lam-1)*(lam+1))
    other = [p-4*d, t, u, 2*u-a, 2*t-b, -d]
    change = dict(zip(coordinates, other))
    multiplier = (phi*sp.Matrix(coordinates)).reshape(3, 3)
    opposite = (phi.subs(lam, -lam)*sp.Matrix(other)).reshape(3, 3)
    assert equal_polynomials(opposite-sp.diag(1, -1, 1)*multiplier)
    assert equal_polynomials(sp.Matrix([4*multiplier.det()-(lam+1)**2*symmetric.det()
                                        +(lam-1)**2*symmetric.xreplace(change).det()]))

    params = [u, -a-b-t, -b-d, a, b, p]
    q = normalized(cusp(True), 1, params, False)
    target = sp.Matrix([[-y, x+y, -x+y], [y+z, -x+y, 0], [-y, x+y+z, -x+z]])
    phi = adjugate_map(q, target, [0, 1, 2, 3, 4, 7], -3)
    assert phi[:, 1].reshape(3, 3).det() == 1

    for aux in (node(0), cusp(False)):
        basis = six_space(aux)
        for m, n in combinations_with_replacement(basis, 2):
            assert bracket(m[1], n[2])+bracket(n[1], m[2]) == sp.zeros(3)
    print('PASS: both singular cubic adjugate identities, nodal pencil, all special parameters and nonfree auxiliaries')


if __name__ == '__main__':
    main()
