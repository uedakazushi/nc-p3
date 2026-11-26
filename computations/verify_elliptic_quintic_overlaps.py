"""Bound the full geometric rank-drop locus of the elliptic quintic cubic overlaps.

All arithmetic before the final characteristic-zero check is in characteristic
181. Norms of interpolated minors bound x over the algebraic closure, not just
at rational points. Local row reduction records the orders used in lifting.
"""
from collections import Counter
from functools import cache
from itertools import combinations, product
from math import factorial

from flint import nmod_mat, nmod_poly
import sympy as sp


PRIME = 181
IMAGINARY = 19
WORDS = list(product(range(4), repeat=3))
POINTS = [None] + [(x, y) for x in range(PRIME) for y in range(PRIME)
                   if (y*y-x**3+x) % PRIME == 0]


@cache
def add(first, second):
    if first is None:
        return second
    if second is None:
        return first
    x, y = first
    a, b = second
    if x == a and (y+b) % PRIME == 0:
        return None
    if first == second:
        slope = (3*x*x-1)*pow(2*y, -1, PRIME) % PRIME
    else:
        slope = (b-y)*pow(a-x, -1, PRIME) % PRIME
    u = (slope*slope-x-a) % PRIME
    return u, (slope*(x-u)-y) % PRIME


def negate(point):
    return None if point is None else (point[0], -point[1] % PRIME)


def multiply(n, point):
    result = None
    for _ in range(n):
        result = add(result, point)
    return result


def project(point):
    if point is None:
        return (0, 0, 0, 1)
    x, y = point
    return (1+IMAGINARY*x*x) % PRIME, x, y, x*y % PRIME


GROUP = [q for q in POINTS if q is None or
         multiply(2, q) == (-q[0] % PRIME, IMAGINARY*q[1] % PRIME)]


def representatives(points):
    result = []
    for q in points:
        if not any(add(q, h) in result for h in GROUP):
            result.append(q)
    return result


def nullvectors(matrix):
    kernel, dimension = matrix.nullspace()
    return [[int(kernel[i, j]) for i in range(matrix.ncols())] for j in range(dimension)]


def character_basis():
    # Recover the projectivity from point pairs and verify it on the whole curve sample.
    rows = []
    for q in POINTS[:12]:
        z = project(q)
        w = project(add(q, GROUP[1]))
        first = next(i for i, value in enumerate(w) if value)
        for i in range(4):
            if i != first:
                rows.append([z[k % 4]*((w[first] if k//4 == i else 0) -
                                      (w[i] if k//4 == first else 0)) % PRIME for k in range(16)])
    entries, = nullvectors(nmod_mat(rows, PRIME))
    matrix = nmod_mat(4, 4, entries, PRIME)
    for q in POINTS:
        z = matrix*nmod_mat(4, 1, project(q), PRIME)
        w = project(add(q, GROUP[1]))
        assert all((int(z[i, 0])*w[j]-int(z[j, 0])*w[i]) % PRIME == 0
                   for i, j in combinations(range(4), 2))
    identity = nmod_mat([[int(i == j) for j in range(4)] for i in range(4)], PRIME)
    eigenvectors = {}
    for n in range(1, PRIME):
        kernel = nullvectors(matrix.transpose()-n*identity)
        if kernel:
            eigenvectors[n], = kernel
    root = next(n for n in range(2, PRIME) if pow(n, 5, PRIME) == 1)
    missing = next(n for n in range(1, PRIME) if n not in eigenvectors and
                   all(n*pow(root, i, PRIME) % PRIME in eigenvectors for i in range(1, 5)))
    return [eigenvectors[missing*pow(root, i, PRIME) % PRIME] for i in range(1, 5)]


BASIS = character_basis()
COSETS = representatives(POINTS)
SAMPLE = (6, 36)


@cache
def coordinates(point):
    return [sum(a*b for a, b in zip(row, project(point))) % PRIME for row in BASIS]


def word_block(weight, delta, shift, fixed):
    words = [w for w in WORDS if (sum(w)+3) % 5 == weight]
    rows, poles = [], []
    for q in fixed:
        u = coordinates(q)
        for right in (False, True):
            offset = add(q, delta) if right else q
            v = coordinates(add(offset, shift))
            for i in range(4):
                rows.append([int((a if right else c) == i) *
                             (u[b]*v[c] if right else u[a]*v[b]) % PRIME for a, b, c in words])
                poles.append(negate(offset))
    return rows, poles


def independent_rows(rows):
    reduced, rank = nmod_mat(rows, PRIME).transpose().rref()
    return [next(i for i in range(reduced.ncols()) if reduced[k, i]) for k in range(rank)]


def monomials(point, degree):
    x, y = point
    return ([pow(x, i, PRIME) for i in range(degree//2+1)] +
            [y*pow(x, i, PRIME) % PRIME for i in range((degree-3)//2+1)])


@cache
def interpolator(degree, omitted):
    points = [q for q in POINTS[1:] if q[0] not in omitted]
    rows = [monomials(q, degree) for q in points]
    selected = independent_rows(rows)
    assert len(selected) == degree
    return [points[i] for i in selected], nmod_mat([rows[i] for i in selected], PRIME).inv()


def minor_norm(weight, delta):
    fixed = COSETS[:3]
    rows, poles = word_block(weight, delta, SAMPLE, fixed)
    selected = independent_rows(rows)
    assert len(selected) == len(rows[0])
    poles = [poles[i] for i in selected]
    finite = Counter(q[0] for q in poles if q is not None)
    degree = 10*sum(finite.values()) + 5*poles.count(None)
    points, inverse = interpolator(degree, tuple(sorted(finite)))

    def cleared_minor(q):
        rows, _ = word_block(weight, delta, q, fixed)
        value = int(nmod_mat([rows[i] for i in selected], PRIME).det())
        for x, multiplicity in finite.items():
            value = value*pow(q[0]-x, 5*multiplicity, PRIME) % PRIME
        return value

    coefficients = inverse*nmod_mat([[cleared_minor(q)] for q in points], PRIME)
    values = [int(coefficients[i, 0]) for i in range(degree)]
    split = degree//2+1
    a, b = nmod_poly(values[:split], PRIME), nmod_poly(values[split:], PRIME)
    for q in POINTS[1:]:
        if q[0] not in finite:
            assert (int(a(q[0]))+q[1]*int(b(q[0]))-cleared_minor(q)) % PRIME == 0
    x = nmod_poly([0, 1], PRIME)
    return a*a-(x**3-x)*b*b


def geometric_candidates(delta):
    norms = [minor_norm(weight, delta) for weight in range(5)]
    bound = norms[0]*norms[1]
    for a, b in combinations(norms, 2):
        bound = bound.gcd(a*b)
    xs = set()
    for factor, _ in bound.factor()[1]:
        assert factor.degree() == 1
        xs.add(-int(factor[0])*pow(int(factor[1]), -1, PRIME) % PRIME)
    assert xs <= {q[0] for q in POINTS[1:]}
    # If clearing a pole cancels its leading term, a determinant vanishing
    # projectively still has zero cleared value there (the pole bound is 5).
    result = []
    for q in [None]+[q for q in POINTS[1:] if q[0] in xs]:
        deficient = []
        for h in range(5):
            rows, _ = word_block(h, delta, q, COSETS[:3])
            if nmod_mat(rows, PRIME).rank() < len(rows[0]):
                deficient.append(h)
        if len(deficient) >= 4:
            result.append(q)
    return result


def coordinate_jets(point):
    x, y = point
    # Derivatives for the invariant vector field Dx=2y, Dy=3x^2-1.
    dx, dy = 2*y, 3*x*x-1
    ddx, ddy = 2*dy, 12*x*y
    columns = [project(point), (2*IMAGINARY*x*dx, dx, dy, dx*y+x*dy),
               (2*IMAGINARY*(dx*dx+x*ddx), ddx, ddy, ddx*y+2*dx*dy+x*ddy)]
    return [[sum(a*b for a, b in zip(row, column))*pow(factorial(i), -1, PRIME) % PRIME
             for i, column in enumerate(columns)] for row in BASIS]


def local_order(weight, delta, point):
    fixed = [q for q in COSETS if add(q, point) is not None and
             add(add(q, delta), point) is not None][:3]
    words = [w for w in WORDS if (sum(w)+3) % 5 == weight]
    rows = []
    for q in fixed:
        u = coordinates(q)
        for right in (False, True):
            v = coordinate_jets(add(point, add(q, delta) if right else q))
            for i in range(4):
                rows.append([[int((a if right else c) == i) *
                              (u[b]*v[c][k] if right else u[a]*v[b][k]) % PRIME
                              for k in range(3)] for a, b, c in words])
    rank = 0
    for col in range(len(words)):
        pivot = next((k for k in range(rank, len(rows)) if rows[k][col][0]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        unit = rows[rank][col]
        for k in range(rank+1, len(rows)):
            numerator = rows[k][col]
            quotient = []
            for h in range(3):
                previous = sum(quotient[g]*unit[h-g] for g in range(h))
                quotient.append((numerator[h]-previous)*pow(unit[0], -1, PRIME) % PRIME)
            for j, old in enumerate(rows[k]):
                basis = rows[rank][j]
                rows[k][j] = [(value-sum(quotient[g]*basis[h-g] for g in range(h+1))) % PRIME
                              for h, value in enumerate(old)]
        rank += 1
    if rank == len(words):
        return 0
    assert rank == len(words)-1
    return min([h for row in rows[rank:] for cell in row for h, v in enumerate(cell) if v]+[3])


def torsion_relations(a):
    x, y = sp.symbols("x y")
    def projective(u, v):
        return [1+sp.I*u*u, u, v, u*v]
    slope = 3*a*a-1
    shifted = projective(a+slope/(x-a), -slope*y/(x-a)**2)
    entries = []
    for u, v in product(projective(x, y), shifted):
        numerator = sp.cancel(u*v*(x-a)**3)
        entries.append(sp.Poly(sp.rem(numerator, y*y-x**3+x, y), x, y))
    monomials = sorted({power for f in entries for power in f.monoms()})
    coefficient_matrix = sp.Matrix([[f.coeff_monomial(power) for f in entries] for power in monomials])
    return sp.Matrix.hstack(*coefficient_matrix.nullspace())


def main():
    assert IMAGINARY**2 % PRIME == PRIME-1 and len(GROUP) == 5
    torsion = [q for q in POINTS if multiply(10, q) is None]
    assert len(torsion) == 100
    deltas = representatives(torsion)
    assert len(deltas) == 20
    for delta in deltas:
        actual = set(geometric_candidates(delta))
        expected = set()
        if delta is None:
            expected = set(GROUP)
        if delta == (0, 0):
            expected = {add(q, h) for q in [(1, 0), (PRIME-1, 0)] for h in GROUP}
        assert actual == expected, (delta, actual, expected)
    for delta, point, expected in [(None, None, 2), ((0, 0), (1, 0), 1),
                                    ((0, 0), (PRIME-1, 0), 1)]:
        assert [local_order(h, delta, point) for h in range(5)] == [0]+[expected]*4
    plus, minus = torsion_relations(1), torsion_relations(-1)
    assert plus.cols == minus.cols == 6
    for a, b in [(plus, minus), (minus, plus)]:
        matrix = sp.kronecker_product(a, sp.eye(4)).row_join(sp.kronecker_product(sp.eye(4), b))
        assert matrix.to_DM().rank() == 44
    print("PASS: 20 torsion offsets, geometric minor bounds, local orders, and both exact torsion overlaps")


if __name__ == "__main__":
    main()
