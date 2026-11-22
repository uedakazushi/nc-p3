"""Degree-one reductions and both golden-ratio boundary disks of the octic at 109."""

from concurrent.futures import ProcessPoolExecutor
from itertools import product
from flint import fq_default_ctx, nmod, nmod_poly
from sympy import sqrt_mod

from rank_fifteen_octic import (
    finite_boundary,
    finite_windows,
    frame_from_traces,
    octic_frame,
)
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile


PRIME = 109
ROOTS = tuple(j for j in range(1, PRIME // 2 + 1) if (j**4 - j * j + 1) % PRIME == 0)


def ordinary_residue(tx):
    t = nmod(tx, PRIME)
    roots = sqrt_mod(int(t**8 + 14 * t**4 + 1), PRIME, all_roots=True)
    if roots:
        assert len(roots) == 2
        ordinates = [nmod(y, PRIME) for y in roots]
        unit = nmod(1, PRIME)

        def wrap(x):
            return x
    else:
        field = fq_default_ctx(PRIME, 2, "a")
        t = field(tx)
        unit = field.one()
        ordinates = [(t**8 + 14 * t**4 + 1).sqrt()]
        wrap = Scalar
    for j, y, color in product(ROOTS, ordinates, range(2)):
        rules = finite_windows(octic_frame(t, y, j * unit, color), unit, wrap)
        assert verify_central_profile(rules, wrap(unit))
    print("Ordinary t orbit", tx, "verified", flush=True)


def trace_boundary(k0):
    field = fq_default_ctx(PRIME, 2, "a")
    k = field(k0)
    assert (k * k + k - 1).is_zero()
    H = (k * k - 1).sqrt()
    v = (4 * k * k - 1).sqrt()
    assert H**PRIME == -H and v**PRIME == -v
    i = field(-1).sqrt()
    assert i**PRIME == i
    for r, j, color in product((i, -i), ROOTS, range(2)):
        frame = frame_from_traces(k, H, r, 2 * k + v, field(j), color)
        rules = finite_windows(frame, field.one(), Scalar)
        assert verify_central_profile(rules, Scalar(field.one()))
    print("Trace boundary K", k0, "verified", flush=True)


def run():
    assert ROOTS == (8, 41)
    bad = {
        int(x)
        for poly in finite_boundary()
        for x, _ in nmod_poly(list(poly), PRIME).roots()
    }
    assert bad == {0, 1, 5, 22, 33, 37, 53, 56, 72, 76, 87, 104, 108}
    valid = set(range(PRIME)) - bad
    representatives = sorted(t for t in valid if t <= -pow(t, -1, PRIME) % PRIME)
    assert len(representatives) == 48
    with ProcessPoolExecutor(max_workers=3) as pool:
        boundary = [pool.submit(trace_boundary, k) for k in (10, 98)]
        for _ in pool.map(ordinary_residue, representatives):
            pass
        for check in boundary:
            check.result()


if __name__ == "__main__":
    run()
