"""Project the opposite determinants from the signed fixed curves over the entire slope base."""

import sympy as sp

from verify_fixed_exchange_isogeny import exchange
from verify_fixed_line_elliptic_curves import bad_sections, ordinary_unit
from verify_fixed_opposite_rank_drop import a, factors, quartics
from verify_klein_elliptic_completion import coeffs, gram_twist
from verify_uniform_fixed_quartic_support import u, z, l
from verify_uniform_fixed_quartic_residuals import companion, r as cover_root


def projected_quartics(sign):
    N, D = exchange(sign, z, l) * sp.Matrix([a, 1])
    A = sp.Matrix([[a, 1], [sign * z, sign * a]])
    e, f = N * u, D * u
    B = sp.Matrix([[e, f], [sign * z * f, sign * e]])
    co, _, _ = coeffs(*A, z, l)
    eq = sp.factor(sign * z * co[1] * f**2 + f - co[0] - sign * co[1] * a * a)
    (curve,) = [g for g in factors(sp.fraction(eq)[0]) if g.has(u)]
    opp = gram_twist(A, B, z, l)[1]
    delta = a * a - z
    excluded = factors(a * delta * N * D)
    projected = []
    for h in range(4):
        indices = [4 * j + (j ^ h) for j in range(4)]
        block = opp.extract(indices, indices).applyfunc(sp.factor)
        det = sp.factor(block.det(method="domain-ge"))
        if h in (2, 3):
            assert sp.cancel(det + (l - 1) ** 2 * (l * z - 1) ** 2 / l**2) == 0
            continue
        if sign == 1:
            # Scale one row or column to expose the two sign representations.
            scaled = block.copy()
            if h == 0:
                scaled[3, :] = scaled[3, :] / z
                m, n = a, e
                unit = -1 / (l * u * z * delta)
                residual = a * a + e * e - 2 * z * f
            else:
                scaled[:, 1] = scaled[:, 1] / z
                m, n = sp.S.One, f
                unit = -1 / (l * u * delta)
                residual = 1 + f * f - 2 * a * e / z
            X, Y = scaled[2, 2], scaled[2, 3]
            model = sp.Matrix([[m, n, 1, 1], [n, m, 1, 1], [1, 1, X, Y], [1, 1, Y, X]])
            assert all(sp.cancel(v) == 0 for v in scaled - model)
            assert sp.cancel(X - Y - unit * (m - n)) == 0
            assert sp.cancel((m + n) * (X + Y) - 4 + unit * residual) == 0
            projection = sp.resultant(curve, m - n, u)
            assert projection and all(
                P in excluded or ordinary_unit(P, z, l) for P in factors(projection)
            )
            target = (
                -((a - e) ** 2)
                * (a * a + e * e - 2 * z * f)
                / (l * l * u * u * z * delta**2)
                if h == 0
                else -z
                * (1 - f) ** 2
                * (1 + f * f - 2 * a * e / z)
                / (l * l * u * u * delta**2)
            )
            assert sp.cancel(det - target) == 0
        num, den = sp.fraction(det)
        assert all(
            g in factors(u * delta) or ordinary_unit(g, z, l) for g in factors(den)
        )
        extra = []
        for g in factors(num):
            projection = sp.resultant(curve, g, u) if g.has(u) else g
            assert projection
            for P, multiplicity in sp.factor_list(projection)[1]:
                if P in excluded or ordinary_unit(P, z, l):
                    continue
                assert sign == 1 and sp.degree(P, a) == 4 and multiplicity == 1
                extra.append(P)
        if sign == 1:
            (P,) = extra
            assert (
                sp.Poly(P.subs({z: 3, l: 2}), a).monic()
                == sp.Poly(quartics[h], a).monic()
            )
            projected.append(P)
        else:
            assert not extra
    if sign == 1:
        for sig, av, ee, ff in bad_sections(z, l):
            if sig == 1:
                assert all(
                    ordinary_unit(value, z, l)
                    for value in (
                        av * av + ee * ee - 2 * z * ff,
                        z * (1 + ff * ff) - 2 * av * ee,
                    )
                )
        for kind in ("zero", "infinity"):
            av, bv, ee, ff = companion(kind)
            for value in (av - ee, bv - ff):
                top, bottom = sp.fraction(sp.factor(value))
                assert ordinary_unit(top.subs(cover_root**2, z), z, l)
                assert bottom != 0
        P, Q = projected
        assert sp.cancel(a**4 * P.subs(a, z / a) / Q) == -(z**2)
        common = sp.factor(sp.resultant(P, Q, a))
        exceptional = [T for T in factors(common) if not ordinary_unit(T, z, l)]
        (T,) = exceptional
        expected = (
            (l * l * z + 1) ** 2
            + 2 * l * (z + 1) * (l * l * z + 1)
            + l * l * (z * z - 14 * z + 1)
        )
        assert sp.expand(T - expected) == 0
    print(
        sign, "opposite determinants and exact projected factors verified", flush=True
    )


if __name__ == "__main__":
    for sign in (1, -1):
        projected_quartics(sign)
