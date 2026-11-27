"""Smooth cubic surfaces and their 27 intersections on the marked 20-point chart."""
from collections import defaultdict
from itertools import combinations, combinations_with_replacement, permutations, product
from math import factorial

import numpy as np

from verify_twenty_mixed_bracket_invariants import (
    DIGITS, MOD, P, configurations, kernel_over_dvr, lift, lift_points, marked_tensor, tensor,
)


def monomials(n, degree):
    return [tuple(indices.count(i) for i in range(n))
            for indices in combinations_with_replacement(range(n), degree)]


def cross(x, y):
    return np.array([x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3]
                     for i in range(3)], dtype=object) % MOD


def determinant(rows):
    return int(rows[0] @ cross(rows[1], rows[2])) % MOD


def minors(rows):
    return np.array([(-1)**j*determinant(np.delete(rows, j, axis=0))
                     for j in range(4)], dtype=object) % MOD


def onto(variables, degree, equations):
    target = {m: i for i, m in enumerate(monomials(variables, degree))}
    columns = []
    for eq, d in equations:
        for h in monomials(variables, degree-d):
            col = [0]*len(target)
            for q, c in eq.items():
                col[target[tuple(i+j for i, j in zip(q, h))]] = c
            columns.append(col)
    return kernel_over_dvr(np.array(columns, dtype=object).T, DIGITS)[1]


def smooth(E):
    cubic = defaultdict(int)
    for perm in permutations(range(3)):
        sign = (-1)**sum(perm[i]>perm[j] for i, j in combinations(range(3), 2))
        for x, y, z in product(range(4), repeat=3):
            keys = (x, y, z)
            cubic[tuple(keys.count(i) for i in range(4))] += (
                sign*E[x, 0, perm[0]]*E[y, 1, perm[1]]*E[z, 2, perm[2]])
    derivatives = []
    for i in range(4):
        partial = {tuple(e-(i==j) for j, e in enumerate(q)): q[i]*c % MOD
                   for q, c in cubic.items() if q[i]}
        derivatives.append((partial, 2))
    return onto(4, 5, derivatives)


def disjoint_base(E, F):
    eqs = []
    for space in (E, F):
        for I in combinations(range(4), 3):
            cubic = defaultdict(int)
            for xyz in product(range(3), repeat=3):
                cubic[tuple(xyz.count(h) for h in range(3))] += determinant(
                    [space[i, j] for i, j in zip(I, xyz)])
            eqs.append((cubic, 3))
    return onto(3, 8, eqs)


def interpolate(values):
    coeff, newton = [0]*len(values), [1]
    for j, y in enumerate(values):
        error = y-sum(c*pow(j, i, MOD) for i, c in enumerate(coeff))
        change = error*pow(factorial(j), -1, MOD) % MOD
        for i, c in enumerate(newton):
            coeff[i] = (coeff[i]+change*c) % MOD
        newton = [(c-j*d) % MOD for c, d in zip([0]+newton, newton+[0])]
    return coeff


def no_common_curve(E, F):
    values = []
    for j in range(11):
        x = np.array([1, j, 2+3*j], dtype=object)
        images = []
        for space in (E, F):
            Z = sum(a*M for a, M in zip(minors(x @ space), space)) % MOD
            images.append(cross(Z[1], Z[2]))
        v = cross(*images)
        values.append([v @ coeff % MOD for coeff in ((1,2,5), (2,3,7))])
    forms = list(map(interpolate, zip(*values)))
    sylvester = [[0]*j+f+[0]*(9-j) for f in forms for j in range(10)]
    return kernel_over_dvr(sylvester, DIGITS)[1]


def surface(t, start):
    axes = tuple((start+i)%4 for i in range(4))
    return t.transpose(axes)[0, 1:, :, 1:].transpose(1, 0, 2)


def vector_and_derivative(E, x, dx, fixed):
    rows, drows = x @ E, dx @ E
    a = minors(rows)
    da = []
    for j in range(4):
        I = [i for i in range(4) if i!=j]
        da.append((-1)**j*sum(determinant([drows[i] if i==s else rows[i] for i in I]) for s in I))
    Z = sum(f*M for f, M in zip(a, E)) % MOD
    dZ = sum(f*M for f, M in zip(da, E)) % MOD
    i, j = (fixed+1)%3, (fixed+2)%3
    return cross(Z[i], Z[j]), (cross(dZ[i], Z[j])+cross(Z[i], dZ[j])) % MOD


def transverse(E, F, marked, points):
    pivot = next(i for i in range(4) if marked[i] % P)
    for b in points:
        vec = b-b[pivot]*pow(int(marked[pivot]), -1, MOD)*marked
        x = np.delete(vec, pivot) % MOD
        if not any(x):
            assert np.array_equal(b, marked)
            continue
        precision = DIGITS
        while not any(x % P):
            x //= P
            precision -= 1
        fixed = next(i for i in range(3) if x[i] % P)
        derivatives = []
        for i, dx in enumerate(np.eye(3, dtype=object)):
            if i == fixed:
                continue
            f, df = vector_and_derivative(E, x, dx, fixed)
            g, dg = vector_and_derivative(F, x, dx, fixed)
            assert any(f % P**precision) and any(g % P**precision)
            assert not any(cross(f, g) % P**precision)
            derivatives.append(cross(df, g)+cross(f, dg))
        assert any(cross(*derivatives) % P**precision)


def check():
    W = tensor(lift())
    odd, first, third = configurations()
    tracks = np.array([lift_points(W, t) for t in zip(first, odd, third, odd)],
                      dtype=object).transpose(1, 0, 2)
    for j in range(20):
        pts = tracks[:, j, :]
        normal = np.vectorize(int, otypes=[object])(marked_tensor(W, pts))
        raw = marked_tensor(W, pts, normalize_gram=False)
        for i in (0, 1):
            E, F = surface(normal, i), surface(normal, i+2).transpose(0, 2, 1)
            assert np.array_equal(E[0], np.eye(3)) and np.array_equal(E[0], F[0])
            assert min(smooth(E), smooth(F), disjoint_base(E,F), no_common_curve(E,F)) > 0
            assert kernel_over_dvr(np.vstack([E, F[1:]]).reshape(7, 9), DIGITS)[1] > 0
            brackets = np.array([(g@f-f@g).ravel()[:8] for g in E[1:] for f in F[1:]], dtype=object)
            assert kernel_over_dvr(brackets.T, DIGITS)[1] > 0
            transverse(surface(raw, i), surface(raw, i+2).transpose(0, 2, 1),
                       pts[i+1], tracks[i+1])
    print('PASS: 80 smooth cubics, 40 proper pairs, transverse marked intersections and bracket ranks')


if __name__ == '__main__':
    check()
