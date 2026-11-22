"""Certify the five-dimensional quartic center on both entire ordinary rank-14 components."""

import sympy as sp
from flint import fmpz_poly, nmod, nmod_poly

from rank_fourteen_function_field import (
    CurveElement,
    excluded_factors,
    irreducible_factors,
    one,
    relations,
)
from verify_fixed_branch_regular import check_window
from verify_klein_center_lift import (
    add_multiple,
    central_equations,
    echelon,
    field_paths,
)
from verify_rank_fourteen_descent import frame, slopes


def path_chart(windows):
    denominators = set()

    def reduce(rows):
        pivots = echelon(rows)
        for row in pivots.values():
            for value in row.values():
                denominators.add(tuple(map(int, value.denominator)))
        return pivots

    basis, image = field_paths(windows, one, reduce)
    extra = (
        set().union(*(irreducible_factors(fmpz_poly(list(d))) for d in denominators))
        - excluded_factors()
    )
    assert sorted(len(p) - 1 for p in extra) == [8]
    return basis, image, extra


def rank_minor(rows, unit=one):
    pivots = {}
    determinant = unit
    for original in rows:
        row = original.copy()
        for j in sorted(pivots):
            if j in row:
                add_multiple(row, pivots[j], -row[j])
        if row:
            j = min(row)
            pivot = row[j]
            determinant *= pivot
            inverse = pivot.inverse()
            pivots[j] = {k: v * inverse for k, v in row.items()}
    return len(pivots), determinant


def certify_residue(polynomial, prime, t0):
    poly = nmod_poly(polynomial, prime)
    assert (
        poly.degree() == len(polynomial) - 1 and not poly(t0) and poly.derivative()(t0)
    )
    t = nmod(t0, prime)
    r, lam = slopes(t)
    z = r * r
    assert (
        z
        * lam
        * (z * z - 1)
        * (lam * lam * z - 1)
        * (lam * lam - 1)
        * (lam * lam * z * z - 1)
    )
    signs = sp.sqrt_mod(-1, prime, all_roots=True)
    roots = sp.sqrt_mod(int(lam * r), prime, all_roots=True)
    assert len(signs) == len(roots) == 2
    for i in signs:
        for s in roots:
            X, b, E, F = frame(nmod(s, prime), nmod(i, prime))
            scale = b if b else X
            aa, bb, ee, ff = map(
                int, (r * X / scale, b / scale, r * E / scale, F / scale)
            )
            zz, ll = sp.Integer(int(z)), sp.Integer(int(lam))
            A = sp.Matrix([[aa, bb], [zz * bb, aa]])
            B = sp.Matrix([[ee, ff], [zz * ff, ee]])
            check_window(
                A,
                B,
                zz,
                ll,
                (len(polynomial) - 1, prime, t0, i, s),
                prime,
                opposite_rank=14,
            )


if __name__ == "__main__":
    basis, image, extra = path_chart(relations())
    expected_degrees = ([8, 8, 20], [8, 34], [8], [8])
    for h, (n, rank) in enumerate(((36, 33), (40, 38), (32, 32), (32, 32))):
        rows, size = central_equations(basis, image, h)
        assert size == n and len(rows) == 224
        computed, det = rank_minor(rows)
        assert computed == rank
        factors = irreducible_factors(det.norm_numerator()) - excluded_factors()
        assert sorted(len(f) - 1 for f in factors) == expected_degrees[h]
        extra.update(factors)
    assert len(extra) == 4
    for poly, (prime, t0) in zip(
        sorted(extra), ((157, 141), (101, 89), (101, 45), (101, 85))
    ):
        certify_residue(list(poly), prime, t0)
    print(
        "Uniform quartic center dimensions (3,2,0,0); all four residual orbits verified.",
        flush=True,
    )
