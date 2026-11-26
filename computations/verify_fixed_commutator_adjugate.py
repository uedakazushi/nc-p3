"""The universal six quadrics of a fixed smooth ordered three-commutator target."""
from itertools import combinations_with_replacement

import sympy as sp
from verify_moore_commutator_doubling import moore
from verify_prescribed_discriminant_seeds import axial, bracket


def sparse_triple(diag, positive, negative):
    out = [sp.diag(*diag), sp.zeros(3), sp.zeros(3)]
    for i in range(3):
        out[1][i, (i+1) % 3] = positive[i]
        out[2][i, (i-1) % 3] = negative[i]
    return out


def rotations(q):
    j = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    return [[j**n*q[(i-n) % 3]*j.T**n for i in range(3)] for n in range(3)]


def linear_basis(a, b, c):
    weights = [c/b, b/a, a/c]
    s = [a/b, c**2/b**2, c/a]
    u = sparse_triple([0]*3, s, [m*h for m, h in zip(weights, s)])
    d = [(2*a**3-b**3-c**3)/(3*a*b**2), (-a**3-b**3+2*c**3)/(3*a*b**2),
         (-a**3+2*b**3-c**3)/(3*a*b**2)]
    s = [-c/a, -a/b, -c**2/b**2]
    w = sparse_triple(d, s, [m*h for m, h in zip(weights, s)])
    return rotations(u)+rotations(w)


def zero(matrix):
    return all(sp.cancel(entry) == 0 for entry in matrix)


def main():
    a, b, c = sp.symbols('a b c', nonzero=True)
    basis = linear_basis(a, b, c)
    v = sp.symbols('x y z')
    matrix = moore((a, b, c), v)
    equations = sp.Matrix.hstack(*(matrix.diff(x) for x in v))
    for q in basis:
        for i, j in combinations_with_replacement(range(3), 2):
            h = sp.zeros(3)
            h[i, j] = h[j, i] = 1
            assert zero(equations*sp.Matrix.vstack(*[axial(h*u-u.T*h) for u in q]))

    doubled = [a*(b**3-c**3), b*(c**3-a**3), c*(a**3-b**3)]
    reference = moore(doubled, v)[[0, 2, 1], :]
    targets = [reference.diff(x) for x in v]
    inv = targets[0].inv()
    coefficients = {}
    for i, j in combinations_with_replacement(range(6), 2):
        cs = []
        for k in range(3):
            n, m = (k+1) % 3, (k+2) % 3
            z = bracket(basis[i][n], basis[j][m])+bracket(basis[j][n], basis[i][m])
            cs.append((z/(1+int(i == j))).applyfunc(sp.factor))
        left = (cs[0]*inv).applyfunc(sp.factor)
        assert all(zero(x-left*y) for x, y in zip(cs[1:], targets[1:]))
        coefficients[i, j] = left

    # Every entry vanishes on Z=p p^t. Group by the 15 monomials in p of degree four.
    p = sp.symbols('p:3')
    z = [x*x for x in p]+[p[(i+1) % 3]*p[(i+2) % 3] for i in range(3)]
    degree_four = {}
    for (i, j), mat in coefficients.items():
        monomial = z[i]*z[j]
        degree_four[monomial] = degree_four.get(monomial, sp.zeros(3))+mat
    assert len(degree_four) == 15 and all(zero(mat) for mat in degree_four.values())

    # The two coefficients from each homogeneous pair have disjoint cyclic diagonals.
    minors = [(a**3-c**3)/(a**4*b**7), -c*(a**3-b**3)/(a**4*b**8),
              -(b**3-c**3)/(a**3*b**8)]
    occupied = set()
    for i in range(3):
        x, y = coefficients[3+i, 3+i], coefficients[i, 3+i]
        ids = [j for j in range(9) if x[j] != 0 or y[j] != 0]
        assert len(ids) == 3 and not occupied.intersection(ids)
        occupied.update(ids)
        j, k = ids[:2]
        assert sp.cancel(x[j]*y[k]-x[k]*y[j]-minors[i]) == 0
    # Recover the square-root Veronese on the left multipliers, with its actual
    # row and column embeddings.  These coefficients are the 00 and 12 cofactors.
    surface = sp.zeros(3)
    for i in range(3):
        surface -= p[i]**2 * coefficients[3+i, 3+i]
        surface -= p[(i+1) % 3]*p[(i+2) % 3]*coefficients[i, 3+i]
    adj = moore((a, b, c), p)[[0, 2, 1], :].adjugate()
    origin = {p[0]: 1, p[1]: 0, p[2]: 0}
    factor = sp.factor(surface[0, 0].subs(origin)/adj[0, 0].subs(origin))
    assert factor != 0 and zero(surface-factor*adj)
    print("PASS: 21 quadrics, Veronese base, multiplier rank six and its Moore adjugate surface")


if __name__ == '__main__':
    main()
