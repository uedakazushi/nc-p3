"""Relative regular charts at the two simple ramification values of the genus-three fiber."""

from concurrent.futures import ProcessPoolExecutor
from itertools import product

from flint import fq_default_ctx, nmod
from verify_rank_fifteen_pencil_ramified_ends import (
    relative_limit,
    trace,
    local_frame,
    vanishes,
    value_grad,
    grassmann_chart,
    minors,
)
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile


def verify(job):
    c0, delta, color = job
    field = fq_default_ctx(109, 2, "a")
    c = field(c0)
    assert (c * c - c + 1).is_zero()
    k, other = trace(c, delta, delta)
    assert k == field(-2 * delta) and not (2 * k + delta).is_zero()
    radical = (other * other - 1).sqrt()
    assert not radical.is_zero()
    for sign, rad in product((-1, 1), (-radical, radical)):
        rules, _ = relative_limit(c, k, other + rad, delta, delta, sign, color, Scalar)
        assert verify_central_profile(rules, Scalar(field.one()))
    print(
        "Special fiber",
        c0,
        "trace",
        -2 * delta,
        "presentation",
        color,
        "verified",
        flush=True,
    )


def verify_node(job):
    """Remove separate r+1 and w-1 factors in the two-dimensional local ring."""
    c0, sign, color = job
    p = 109
    c = nmod(c0, p)
    one = nmod(1, p)
    k = nmod(0, p)
    r = -one
    w = one
    n = sign * c
    assert c * c - c + 1 == 0 and n == sign * (c - 1) / c
    alternatives = []
    for delta in (1, -1):
        eps = -delta
        x = w if delta == 1 else r
        alternatives.append(local_frame(c, k, n, eps * one, x, delta, color))
    frame = alternatives[0][0]
    assert frame == alternatives[1][0]
    ws = []
    counts = [0, 0]
    for side in minors():
        rules = []
        for h, group in side:
            evaluations = [
                {cols: value_grad(g, frame, vel) for cols, g in group.items()}
                for _, vel in alternatives
            ]
            L = {cols: v for cols, (v, _) in evaluations[0].items() if v}
            if not L:
                for idx, delta in enumerate((1, -1)):
                    if all(
                        vanishes(g, delta, -delta, sign, color) for g in group.values()
                    ):
                        L = {cols: v for cols, (_, v) in evaluations[idx].items() if v}
                        if L:
                            counts[idx] += 1
                            break
            assert L, (job, len(ws), h)
            for row in grassmann_chart(L, one):
                rules.append((h, {(a, a ^ h): v for a, v in row.items()}))
        ws.append(rules)
    assert counts == [2, 2]
    assert verify_central_profile(ws, one)
    print("NODE", job, counts, flush=True)


if __name__ == "__main__":
    roots = [c for c in range(109) if (c * c - c + 1) % 109 == 0]
    assert len(roots) == 2
    with ProcessPoolExecutor(max_workers=2) as pool:
        list(pool.map(verify, product(roots, (-1, 1), range(2))))
        list(pool.map(verify_node, product(roots, (-1, 1), range(2))))
