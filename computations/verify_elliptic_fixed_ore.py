"""Elliptic smooth-base Ore corrections for automorphisms with fixed points."""

import sympy as sp
from flint import nmod_mat


def add(p, q, prime):
    if p is None:
        return q
    if q is None:
        return p
    x, y = p
    v, w = q
    if x == v:
        if (y + w) % prime == 0:
            return None
        m = 3 * x * x * pow(2 * y, -1, prime) % prime
    else:
        m = (w - y) * pow(v - x, -1, prime) % prime
    z = (m * m - v - x) % prime
    return z, (m * (x - z) - y) % prime


def points(prime):
    return [
        (x, y)
        for x in range(prime)
        for y in range(prime)
        if (y * y - x**3 - 1) % prime == 0
    ]


def kernel(mat):
    m, r = mat.nullspace()
    return [[int(m[i, j]) for i in range(m.nrows())] for j in range(r)]


def rank(prime, points, p, omega):
    rows = []
    for pt in [None, *points]:
        q = add(p, None if pt is None else (omega * pt[0] % prime, pt[1]), prime)
        left = (0, 0, 1) if pt is None else (1, *pt)
        right = (0, 0, 1) if q is None else (1, *q)
        rows.append([a * b % prime for a in left for b in right])
    s = kernel(nmod_mat(rows, prime))
    assert len(s) == 3
    relations = sp.Matrix(s).T
    pure = sp.kronecker_product(relations, sp.eye(3)).row_join(
        sp.kronecker_product(sp.eye(3), relations)
    )
    normals = sp.Matrix(
        kernel(nmod_mat(pure.cols, pure.rows, [int(v) % prime for v in pure.T], prime))
    )
    assert normals.shape == (10, 27)
    echelon, _ = nmod_mat(s, prime).rref()
    pivots = {next(j for j in range(9) if echelon[i, j]) for i in range(3)}
    columns = []
    for i in range(3):
        for j in range(9):
            if j in pivots:
                continue
            op = sp.eye(9)[:, j] * sp.eye(3)[i, :]
            col = (
                normals
                * (
                    sp.kronecker_product(op, sp.eye(3))
                    + omega * sp.kronecker_product(sp.eye(3), op)
                )
                * relations
            )
            columns.append([int(v) % prime for v in col])
    return nmod_mat(columns, prime).rank()


def inversion():
    c = sp.symbols("c", nonzero=True)
    eye = sp.eye(9)
    swap = sp.kronecker_product(sp.eye(3), sp.eye(3)[:, [1, 0, 2]])
    rels = swap * sp.Matrix.hstack(
        *[
            c * eye[:, 4 * j]
            + eye[:, 3 * ((j + 1) % 3) + (j + 2) % 3]
            + eye[:, 3 * ((j + 2) % 3) + (j + 1) % 3]
            for j in range(3)
        ]
    )
    reps = swap * eye[:, [1, 2, 3, 5, 6, 7]]
    cubic = sp.kronecker_product(rels, sp.eye(3)).row_join(
        sp.kronecker_product(sp.eye(3), rels)
    )
    from sympy.polys.matrices import DomainMatrix

    normals = DomainMatrix.from_Matrix(cubic.T).to_field().nullspace().to_Matrix()
    assert sp.lcm_list([sp.factor(a).as_numer_denom()[1] for a in normals]) == 1
    columns = []
    for i in range(3):
        for j in range(6):
            op = reps[:, j] * sp.eye(3)[i, :]
            cols = (
                normals
                * (
                    sp.kronecker_product(op, sp.eye(3))
                    - sp.kronecker_product(sp.eye(3), op)
                )
                * rels
            )
            columns.append(cols.reshape(30, 1))
    pencil = sp.Matrix.hstack(*columns).applyfunc(sp.factor)
    assert DomainMatrix.from_Matrix(pencil).rank() == 13
    rows = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 21, 24]
    cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 15]
    det = pencil.extract(rows, cols).det(method="domain-ge")
    assert sp.factor(det) == 2 * c**60 * (c - 1) ** 28 * (c**2 + c + 1) ** 28
    print("PASS: entire smooth inversion family, rank 13")


def fourth_order():
    from ternary_ore_parameters import presentation_pencil

    x, y = sp.symbols("x y")
    basis = [1, x, y]
    target = [x * x, x, sp.I * y]
    pairs = [
        sp.Poly(sp.rem(a * b, y * y - x**3 + x, y), x, y) for a in basis for b in target
    ]
    monomials = sorted(set().union(*(p.monoms() for p in pairs)))
    rels = sp.Matrix.hstack(
        *sp.Matrix(
            [[p.coeff_monomial(m) for p in pairs] for m in monomials]
        ).nullspace()
    )
    left, right, _, _ = presentation_pencil(rels)
    # dx/y has multiplier i under (x,y)->(-x,i*y); translation does not change it.
    assert (left + sp.I * right).rank() == 15
    print("PASS: fourth-order linear part and its conjugate, zero corrections")


def main():
    inversion()
    fourth_order()
    prime = 109
    curve = points(prime)
    candidates = [
        pt
        for pt in curve
        if add(pt, add(pt, pt, prime), prime) in ((0, 1), (0, prime - 1))
    ]
    assert len(candidates) == 18 and len(curve) == 107
    for omega in [45, 63]:
        assert (omega**2 + omega + 1) % prime == 0
        for pt in candidates:
            assert rank(prime, curve, pt, omega) == 15
    print("PASS: both order-three linear parts, all 18 translations modulo 109")


if __name__ == "__main__":
    main()
