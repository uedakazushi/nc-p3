"""Saturate the quintic rank ideals on the two coordinate fibers and two companion families."""

from functools import cache
import random

import sympy as sp
from flint import nmod_poly

from verify_fixed_geometric_quartic_locus import a, e, f, minor_indices, symbolic_blocks
from verify_fixed_quintic_locus import middle_relations
from verify_uniform_fixed_quartic_residuals import (
    b,
    r,
    companion,
    ordinary_factors,
    scalar_factors,
)
from verify_uniform_fixed_quartic_support import z, l, polynomial, cleared_rows
from verify_varying_fixed_quintics import determinant_factors


def factors_on_family(rows, matrix, allowed, kind, start, character):
    @cache
    def snapshot(rr, lam):
        return [[int(v(0, 0, rr, lam)) % 307 for v in row] for row in rows]

    def factors(rr, lam, seed):
        I, J = minor_indices(matrix, snapshot(rr, lam), seed)
        if len(I) != 86:
            return None
        assert len(J) == 86
        det = determinant_factors([[rows[i][j] for j in J] for i in I])
        return {str(g): g for g, _ in det if str(g) not in allowed}

    first = factors(2, 3, 1)
    assert first is not None
    denominator_factors = [
        polynomial(sp.sympify(key.replace("^", "**"), locals={"z": z, "l": l}))
        for key in allowed
    ]
    for name, P in first.items():
        possible = None
        used = 0
        trials = list(range(2, 307))
        random.Random(2).shuffle(trials)
        for rr in trials:
            terms = [0] * (P.degrees()[3] + 1)
            for (_, _, j, k), c in P.to_dict().items():
                terms[k] = (terms[k] + int(c) * pow(rr, j, 307)) % 307
            poly = nmod_poly(terms, 307)
            if not poly:
                continue
            for lam, _ in poly.roots():
                lam = int(lam)
                if any(int(g(0, 0, rr, lam)) % 307 == 0 for g in denominator_factors):
                    continue
                other = factors(rr, lam, 1 + used % 3)
                if other is None:
                    continue
                assert name not in other
                used += 1
                vertical = set()
                for Q in other.values():
                    res = P.resultant(Q, 3) if Q.degrees()[3] else Q
                    assert res
                    vertical.update(str(g) for g, _ in res.factor()[1])
                vertical -= allowed
                possible = vertical if possible is None else possible & vertical
                if not possible or used >= 5:
                    break
            if possible == set() or used >= 5:
                break
        assert possible == set(), (kind, start, character, name, possible)


def check(kind, blocks):
    av, bv, ee, ff = companion(kind)
    det_a = av * av - r * r * bv * bv
    det_b = ee * ee - r * r * ff * ff
    allowed = ordinary_factors()
    for v in (
        det_a,
        det_b,
        det_a + det_b,
        av * av - r * r * ff * ff,
        ee * ee - r * r * bv * bv,
    ):
        allowed |= scalar_factors(v)
    subs = {a: av, b: bv, e: ee, f: ff, z: r * r}
    for start in (0, 1):
        for h in (0, 2):
            matrix = middle_relations(blocks, start, h)
            values = {
                v: sp.factor(v.subs(subs, simultaneous=True)).subs(r, z)
                for v in set(matrix)
            }
            matrix = matrix.applyfunc(values.get)
            rows = cleared_rows(
                matrix, lambda d: scalar_factors(d.subs(z, r)) <= allowed
            )
            factors_on_family(rows, matrix, allowed, kind, start, h)
        print(kind, start, "no admissible common residual quintic zero", flush=True)


if __name__ == "__main__":
    blocks = symbolic_blocks(1, z, l, b)
    for kind in ("zero", "infinity", 1, -1):
        check(kind, blocks)
