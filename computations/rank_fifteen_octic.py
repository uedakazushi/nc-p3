"""Frame and ordinary divisor on the genus-three rank-fifteen pencil levels."""

from rank_fourteen_function_field import (
    VARIABLE,
    evaluate_relations,
    irreducible_factors,
)
from verify_fixed_geometric_quartic_locus import a, e, f
from verify_uniform_fixed_quartic_residuals import b
from verify_uniform_fixed_quartic_support import z, l


def frame_from_traces(k, H, r, w, j, color):
    zz = r * r
    n = H / j
    aa, bb, ee, ff = 1 + n, (k + j * H) / r, 1 - n, (k - j * H) / r
    if color:
        aa, bb, ee, ff = zz * bb, aa, zz * ff, ee
    return zz, w / r, aa, bb, ee, ff


def octic_frame(t, y, j, color):
    d = 2 * t * (t * t - 1)
    k = -((t * t + 1) ** 2) / d
    H = y / d
    A = k * (t * t - 2 * t - 1) / (t * t + 1)
    B = k * (1 - 2 * t - t * t) / (t * t + 1)
    return frame_from_traces(k, H, k + H * H + H * A, k - H * H + H * B, j, color)


def finite_windows(frame, unit, wrap):
    rules = evaluate_relations(
        dict(zip((z, l, a, b, e, f), frame)), unit, lambda n: n * unit
    )
    return [
        [(h, {word: wrap(c) for word, c in row.items()}) for h, row in side]
        for side in rules
    ]


def finite_boundary():
    t = VARIABLE
    d = 2 * t * (t * t - 1)
    numerator = -((t * t + 1) ** 2)
    factors = [d, numerator**2 + numerator * d - d * d]
    factors.extend(numerator - k * d for k in range(-2, 3))
    return set().union(*(irreducible_factors(g) for g in factors))
