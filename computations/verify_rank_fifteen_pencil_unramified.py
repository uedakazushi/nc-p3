"""Canonical shift quotient and 24 unramified finite ends of a general genus-five level."""

from concurrent.futures import ProcessPoolExecutor
from itertools import product

from flint import nmod
from sympy import symbols, sqrt_mod

from rank_fifteen_octic import finite_windows
from verify_rank_fourteen_geometric_reduction import verify_central_profile


def verify_plane_quartic():
    k, t, s = symbols("K T s")
    h2 = k**2 - t**2
    a2, b2 = ((k + sign * s * t) ** 2 - t**2 for sign in (1, -1))
    d = a2 + b2 - 2 * h2
    assert (s**2 * (a2 - b2) ** 2 - 4 * d * (d + 2 * s**2 * h2)).expand() == 0


def verify_boundary(job):
    zero_trace, color = job
    prime = 101 if zero_trace else 197
    one = nmod(1, prime)
    c, length = 3 * one / 4, 5 * one / 4
    s = c + 1 / c
    assert length**2 == 1 + c * c

    def roots(value):
        values = [nmod(v, prime) for v in sqrt_mod(int(value), prime, all_roots=True)]
        assert len(values) == 2 and values[0] != values[1]
        return values

    traces = [0 * one] if zero_trace else [(-s + v) / 2 for v in roots(s * s + 4)]
    count = 0
    for k in traces:
        r_trace, w_trace = k + (k * k - 1) / s, k - (k * k - 1) / s
        assert k == 0 if zero_trace else r_trace == 0
        for height, dr, dw in product(
            roots(k * k - 1), roots(r_trace**2 - 1), roots(w_trace**2 - 1)
        ):
            r, w = r_trace + dr, w_trace + dw
            z, n = r * r, height / length
            assert z * w * (z - 1) * (w * w - 1)
            a, b, e, f = 1 + n, (k + c * n) / r, 1 - n, (k - c * n) / r
            if color:
                a, b, e, f = z * b, a, z * f, e
            rules = finite_windows((z, w / r, a, b, e, f), one, lambda x: x)
            assert verify_central_profile(rules, one)
            count += 1
    assert count == (8 if zero_trace else 16)
    print(
        "Prime",
        prime,
        "presentation",
        color,
        "trace residues",
        list(map(int, traces)),
        "verified",
        count,
        flush=True,
    )


if __name__ == "__main__":
    verify_plane_quartic()
    with ProcessPoolExecutor(max_workers=2) as pool:
        list(pool.map(verify_boundary, product((True, False), range(2))))
