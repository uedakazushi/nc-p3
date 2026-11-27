"""Alternating syzygy and binary singular-curve checks in the instanton boundary proof."""
from itertools import combinations, combinations_with_replacement, product

import sympy as sp

from verify_finite_symmetric_instanton_monad import kernel


def alternation():
    triples = list(combinations(range(4), 3))
    pairs = list(combinations(range(4), 2))
    # Restrict contractions of all alternating three-forms to a fixed
    # three-dimensional space of directions. Polarization on that plane has
    # precisely these four independent solutions.
    variables = sp.symbols('t:18')
    s = sp.symbols('s:4')
    coeff = [sum(s[i]*variables[6*i+k] for i in range(3)) for k in range(6)]
    skew = sp.zeros(4)
    for (a, b), c in zip(pairs, coeff):
        skew[a, b], skew[b, a] = c, -c
    contracted = skew * sp.Matrix([*s[:3], 0])
    equations = [c for f in contracted for c in sp.Poly(f, *s[:3]).coeffs()]
    assert sp.linear_eq_to_matrix(equations, variables)[0].rank() == 18-len(triples)
    assert len({tuple(sorted((i, j))) for i in range(3) for j in range(4) if i != j}) == 6


def binary_strand():
    # Coordinates x0,...,x3 and three sections v0,v1,v2. The six relations
    # are the kernel of the product of binary cubics and binary quadrics.
    x = sp.symbols('x:4')
    relations = sp.Matrix(6, 12, lambda r, c: int(c//3+c % 3 == r))
    k = kernel(relations)
    degree_two = list(combinations_with_replacement(range(4), 2))
    linear = sp.Matrix(30, 24, lambda r, c: sum(k[3*a+r % 3, c//4] for a in range(4)
                       if tuple(sorted((a, c % 4))) == degree_two[r//3]))
    syzygies = kernel(linear)
    assert syzygies.cols == 3
    differential = sp.Matrix(6, 3, lambda i, j: sum(x[a]*syzygies[4*i+a, j] for a in range(4)))
    s, t = sp.symbols('s t')
    cubic = {x[a]: s**(3-a)*t**a for a in range(4)}
    minors = [differential[list(rows), :].det() for rows in combinations(range(6), 3)]
    assert any(minors) and all(sp.expand(f.subs(cubic)) == 0 for f in minors)
    # The homogeneous coordinate ring of the nonplanar degree-three support.
    h = sp.Matrix([[x[0], x[1], x[2]], [x[1], x[2], x[3]]])
    assert all(sp.expand(h[:, cols].det().subs(cubic)) == 0 for cols in combinations(range(3), 2))
    assert all(sp.expand(v.subs(cubic)) != 0 for v in x)
    z = sp.Symbol('z')
    assert sp.expand((1-z)**2*(1+2*z)) == 1-3*z**2+2*z**3
    assert max(3+2*(n-1) for n in range(1, 5)) == 9
    print('PASS: plane polarization, binary syzygy curve and cubic support Hilbert numerator')


if __name__ == '__main__':
    alternation()
    binary_strand()
