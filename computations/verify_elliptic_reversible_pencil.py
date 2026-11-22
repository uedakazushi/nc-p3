#!/usr/bin/env python3
"""A rational marked elliptic pencil from eight points, residue pairing and Gale duality."""
from itertools import combinations_with_replacement, product

import sympy as sp


def quadrics(points):
    columns = [points[:, i].multiply_elementwise(points[:, j])
               for i, j in combinations_with_replacement(range(4), 2)]
    evaluation = sp.Matrix.hstack(*columns)
    pencil = []
    for coeffs in evaluation.nullspace():
        q = sp.zeros(4)
        for (i, j), entry in zip(combinations_with_replacement(range(4), 2), coeffs):
            q[i, j] = q[j, i] = entry if i == j else entry/2
        pencil.append(q)
    assert len(pencil) == 2
    t = sp.symbols("t")
    discriminant = sp.Poly((pencil[0]+t*pencil[1]).det(method="domain-ge"), t)
    # A reduced binary quartic is the smooth intersection criterion for a quadric pencil.
    assert discriminant.degree() >= 3 and sp.gcd(discriminant, discriminant.diff()).degree() == 0
    coeffs = [discriminant.nth(i) for i in reversed(range(5))]
    a, b, c, d, e = coeffs
    i = 12*a*e-3*b*d+c**2
    j = 72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
    assert 4*i**3 != j**2
    return evaluation, j**2/i**3


def verify_singular_pencil(nodal):
    # Degree-n polynomials on the normalization. Nodal gluing identifies the
    # constant and leading coefficients. Cuspidal gluing prescribes the first jet.
    t = sp.symbols("t")
    def sections(n, parameter):
        if nodal:
            return [parameter+t**n, *[t**k for k in range(1, n)]]
        return [1+parameter*t, *[t**k for k in range(2, n+1)]]
    def coefficients(polynomials, n):
        polynomials = [sp.Poly(g, t) for g in polynomials]
        return sp.Matrix(n+1, len(polynomials), lambda i, j: polynomials[j].nth(i))
    left_parameter, right_parameter = (2, 3) if nodal else (1, 2)
    n_parameter = 6 if nodal else 3
    left, right = sections(4, left_parameter), sections(4, right_parameter)
    middle = sections(8, n_parameter)
    pencil = [middle[0] if nodal else middle[0]+t**8, middle[1]]
    multiples = coefficients([w*g for w in pencil for g in middle], 16)
    assert multiples.rank() == 15
    trace = next(v for v in multiples.T.nullspace()
                 if (v.T*coefficients([a*b for a in middle for b in middle], 16)).rank())
    def pairing(first, second):
        return sp.Matrix(len(first), len(second),
                         lambda i, j: (trace.T*coefficients([first[i]*second[j]], 16))[0])
    mixed = [a*b for a in left for b in right]
    left_square = [a*b for i, a in enumerate(left) for b in left[i:]]
    right_square = [a*b for i, a in enumerate(right) for b in right[i:]]
    assert pairing(mixed, mixed).rank() == 6
    assert pairing(left_square, right_square).rank() == 8


def verify_cycle_degree_obstruction():
    for n in (2, 3, 4):
        degrees = [d for d in product(range(1, 4), repeat=n) if sum(d) == 4]
        for d, e in product(degrees, repeat=2):
            witnesses = []
            for start in range(n):
                for length in range(1, n):
                    arc = {(start+i) % n for i in range(length)}
                    left = sum(d[i] for i in arc) - 1
                    right = sum(e[i] for i in range(n) if i not in arc) - 1
                    witnesses.append(left*right)
            assert (max(witnesses) > 1) == (d != e)


def main():
    x = sp.eye(4).col_join(sp.Matrix([[1,2,3,-6], [2,-1,4,-5], [3,5,-2,-6], [5,1,2,-8]]))
    y = sp.Matrix.hstack(*x.T.nullspace())
    mixed = sp.Matrix.hstack(*(x[:, i].multiply_elementwise(y[:, j])
                               for i in range(4) for j in range(4)))
    assert x.T*y == sp.zeros(4) and mixed.rank() == 7
    xx, invariant_x = quadrics(x)
    yy, invariant_y = quadrics(y)
    assert invariant_x == invariant_y
    section = sp.Matrix([1, 2, 3, 4, 5, 6, 7, -28])
    assert sum(section) == 0 and all(section)
    trace = sp.diag(*(1/v for v in section))
    assert (mixed.T*trace*mixed).rank() == 6
    assert (xx.T*trace*yy).rank() == 8
    e = sp.ones(4, 1)
    f, = y[:4, :].nullspace()
    assert (x*e).multiply_elementwise(y*f) == sp.zeros(8, 1)
    verify_cycle_degree_obstruction()
    for nodal in (True, False):
        verify_singular_pencil(nodal)
    print("PASS nodal and cuspidal normalization controls; smooth apolar quartics with equal discriminant invariant, marked zero, adjacent rank 6 and opposite rank 8")


if __name__ == "__main__":
    main()
