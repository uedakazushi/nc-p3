"""Regular weighted limits at the eight infinite ends of the rank-drop pencil."""

from functools import cache
from itertools import product

import sympy as sp
from flint import nmod
from verify_fixed_geometric_quartic_locus import a, e, f, symbolic_blocks
from verify_rank_fifteen_boundary_disks import WEIGHTS
from verify_rank_fourteen_geometric_reduction import verify_central_profile
from verify_uniform_fixed_quartic_residuals import b
from verify_uniform_fixed_quartic_support import l, z

v = sp.Symbol("v")


@cache
def weighted_window(color, large_r, large_w):
    c, L = (v - 1 / v) / 2, (v + 1 / v) / 2
    scale = 2 * c / L**2
    dr, R = (-2, scale) if large_r else (2, 1 / scale)
    dw, W = (-2, -scale) if large_w else (2, -1 / scale)
    # (order in 1/K, leading coefficient) for a,b,e,f,z,lambda.
    frame = [
        (-1, 1 / L),
        (-1 - dr, (1 + c / L) / R),
        (-1, -1 / L),
        (-1 - dr, (1 - c / L) / R),
        (2 * dr, R * R),
        (dw - dr, W / R),
    ]
    if color:
        aa, bb, ee, ff, zz, ll = frame
        frame = [
            (zz[0] + bb[0], zz[1] * bb[1]),
            aa,
            (zz[0] + ff[0], zz[1] * ff[1]),
            ee,
            zz,
            ll,
        ]
    units = {v, c, L, *[val for _, val in frame]}

    @cache
    def initial(poly):
        coefficients = {}
        for exponents, coefficient in sp.Poly(poly, a, b, e, f, z, l).terms():
            order = sum(degree * k for degree, (k, _) in zip(exponents, frame))
            term = coefficient * sp.prod(
                x**degree for degree, (_, x) in zip(exponents, frame)
            )
            coefficients[order] = coefficients.get(order, 0) + term
        order = min(coefficients)
        value = sp.factor(coefficients[order])
        assert value != 0
        units.add(value)
        return order, value

    weights = WEIGHTS["pole-" if large_r else "pole+", color, 1 if large_w else -1]
    out = []
    for p, side in enumerate(symbolic_blocks(1, z, l, b)):
        rules = []
        for h, (rows, _) in enumerate(side):
            for row in rows:
                terms = {}
                for j, coefficient in enumerate(row):
                    if coefficient == 0:
                        continue
                    (vn, cn), (vd, cd) = map(
                        initial, sp.fraction(sp.cancel(coefficient))
                    )
                    word = j, j ^ h
                    weight = weights[p][j] + weights[(p + 1) % 4][j ^ h]
                    terms[word] = (2 * (vn - vd) + weight, cn / cd)
                minimum = min(order for order, _ in terms.values())
                leading = {
                    word: c for word, (order, c) in terms.items() if order == minimum
                }
                first = next(iter(leading.values()))
                units.add(sp.factor(first))
                rules.append(
                    (h, {word: sp.factor(k / first) for word, k in leading.items()})
                )
        out.append(rules)
    return out, units


def verify(color, large_r, large_w, prime, value):
    rules, units = weighted_window(color, large_r, large_w)

    def residue(expression):
        scalar = expression.subs(v, value)
        assert scalar.is_Rational and int(scalar.q) % prime
        return nmod(int(scalar.p) % prime, prime) / int(scalar.q % prime)

    assert all(residue(unit) for unit in units)
    specialized = [
        [(h, {word: residue(k) for word, k in row.items()}) for h, row in side]
        for side in rules
    ]
    assert verify_central_profile(specialized, nmod(1, prime))


if __name__ == "__main__":
    roots = [j for j in range(109) if (j**4 - j * j + 1) % 109 == 0]
    assert len(roots) == 4
    for color, large_r, large_w in product(range(2), repeat=3):
        verify(color, large_r, large_w, 101, 2)
        for j in roots:
            verify(color, large_r, large_w, 109, (j * j + j) % 109)
        print(color, large_r, large_w, "five split regular limits verified", flush=True)
