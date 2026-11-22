"""Sixteen ramified ends with a single ramified slope on a general genus-five level."""

from itertools import product, combinations
from functools import cache
from math import prod
import sympy as sp
from flint import nmod
from flint_rational_function import rational_function_field
from verify_klein_center_lift import field_paths
from rank_fourteen_function_field import positive_relation_blocks
from verify_rank_fifteen_octic_boundary import grassmann_chart
from verify_rank_fourteen_geometric_reduction import verify_central_profile
from verify_fixed_geometric_quartic_locus import a, e, f
from verify_uniform_fixed_quartic_residuals import b
from verify_uniform_fixed_quartic_support import z, l
from concurrent.futures import ProcessPoolExecutor


@cache
def minors():
    """One primitive polynomial tuple of Pluecker coordinates per character block."""
    out = []
    for side in positive_relation_blocks():
        blocks = []
        for h, (rows, _) in enumerate(side):
            m = sp.Matrix(rows)
            values = {
                cols: sp.cancel(m[:, list(cols)].det())
                for cols in combinations(range(4), len(rows))
            }
            den = sp.lcm([sp.denom(v) for v in values.values()])
            group = {cols: sp.cancel(v * den) for cols, v in values.items() if v}
            common = sp.gcd(list(group.values()))
            group = {
                cols: sp.Poly(sp.cancel(v / common), z, l, a, b, e, f)
                for cols, v in group.items()
            }
            blocks.append((h, group))
        out.append(blocks)
    return out


def trace(c, delta, eps):
    s = c + 1 / c
    k = -delta * s - eps
    other = 2 * k - eps
    return k, other


def phases(c, k, n, r, w, color):
    zz = r * r
    l = w / r
    a, b, e, f = 1 + n, (k + c * n) / r, 1 - n, (k - c * n) / r
    if color:
        a, b, e, f = zz * b, a, zz * f, e
    return zz, l, a, b, e, f


def local_frame(p, c, k, n, eps, x, delta, color):
    unit = nmod(1, p)
    zero = nmod(0, p)
    r, w = (eps, x) if delta == 1 else (x, eps)
    dr, dw = (unit, zero) if delta == 1 else (zero, unit)
    zz, ll, a, b, e, f = phases(c, k, n, r, w, 0)
    vz, vl, va, vb, ve, vf = (
        2 * r * dr,
        dw / r - w * dr / (r * r),
        zero,
        -b * dr / r,
        zero,
        -f * dr / r,
    )
    if color:
        vz, vl, va, vb, ve, vf = vz, vl, vz * b + zz * vb, va, vz * f + zz * vf, ve
    return phases(c, k, n, r, w, color), (vz, vl, va, vb, ve, vf)


def value_grad(poly, point, tangent):
    val = grad = 0
    for powers, coeff in poly.terms():
        val += int(coeff) * prod(x**n for x, n in zip(point, powers))
        for i, (n, d) in enumerate(zip(powers, tangent)):
            if n and d:
                grad += (
                    int(coeff)
                    * n
                    * d
                    * prod(
                        x ** (m - (j == i))
                        for j, (x, m) in enumerate(zip(point, powers))
                    )
                )
    return val, grad


C, Q = sp.symbols("c q")


@cache
def vanishes(poly, delta, eps, sign, color):
    """Test the constant coefficient identically on the quadratic boundary divisor."""
    s = C + 1 / C
    k, other = trace(C, delta, eps)
    n = sign * (C + delta * eps) / C
    r, w = (eps, Q) if delta == 1 else (Q, eps)
    frame = phases(C, k, n, r, w, color)
    # In the quadratic boundary algebra Q^2-2*other*Q+1=0, Q is a unit.
    numerator, denominator = sp.fraction(
        sp.cancel(
            poly.as_expr().subs(dict(zip((z, l, a, b, e, f), frame)), simultaneous=True)
        )
    )
    modulus = sp.Poly(Q * Q - 2 * other * Q + 1, Q, domain=sp.QQ.frac_field(C))
    assert sp.Poly(denominator, Q, domain=modulus.domain).gcd(modulus).degree() == 0
    return sp.Poly(numerator, Q, domain=modulus.domain).rem(modulus).is_zero


def verify(job):
    delta, eps, color = job
    p = 101
    c = nmod(3, p) / 4
    k, T = trace(c, delta, eps)
    root = sp.sqrt_mod(int(T * T - 1), p, all_roots=True)
    assert len(root) == 2 and 2 * k + (c + 1 / c) * delta
    count = 0
    jumps = 0
    for sign, rad in product((-1, 1), root):
        n = sign * (c + delta * eps) / c
        x = T + nmod(rad, p)
        frame, vel = local_frame(p, c, k, n, nmod(eps, p), x, delta, color)
        ws = []
        for side in minors():
            rules = []
            for h, group in side:
                pair = {cols: value_grad(g, frame, vel) for cols, g in group.items()}
                L = {cols: x for cols, (x, y) in pair.items() if x}
                if not L:
                    assert all(
                        vanishes(poly, delta, eps, sign, color)
                        for poly in group.values()
                    )
                    L = {cols: y for cols, (x, y) in pair.items() if y}
                    jumps += 1
                assert L
                for row in grassmann_chart(L, nmod(1, p)):
                    rules.append((h, {(a, a ^ h): v for a, v in row.items()}))
            ws.append(rules)
        assert verify_central_profile(ws, nmod(1, p))
        count += 1
    assert count == 4 and jumps == 8
    print(
        "Branch",
        delta,
        eps,
        "color",
        color,
        "trace residue",
        int(k),
        "verified",
        count,
        flush=True,
    )


@cache
def first(poly):
    cs = poly.terms()
    return sum(v for m, v in cs), tuple(sum(m[i] * v for m, v in cs) for i in range(6))


def verify_trace_jump(job):
    field = rational_function_field(("c",))
    k, sign1, sign2, color = job
    dr = sign1 * (C + k)
    dw = sign2 * (C - k)
    dz = 2 * k * dr
    dl = (dw - dr) * k
    db = (C - dr) * k
    df = (-C - dr) * k
    tangent = (
        (dz, dl, 1, db, -1, df) if color == 0 else (dz, dl, dz + db, 1, dz + df, -1)
    )
    ws = []
    flags = []
    for p, side in enumerate(minors()):
        rules = []
        for h, cs in side:
            values = {J: first(poly) for J, poly in cs.items()}
            group = {J: field.from_sympy(x) for J, (x, deriv) in values.items() if x}
            if not group:
                flags.append((p, h))
                group = {
                    J: v
                    for J, (_, d) in values.items()
                    if (v := field.from_sympy(sum(a * b for a, b in zip(d, tangent))))
                }
            assert group
            for row in grassmann_chart(group, field.one):
                rules.append((h, {(a, a ^ h): v for a, v in row.items()}))
        ws.append(rules)
    bs, _ = field_paths(ws, field.one)
    dims = [
        [sum(len(bs[p, n, h]) for h in range(4)) for p in range(4)] for n in range(3, 6)
    ]
    depends = any(
        v.n.total_degree() > 0 or v.d.total_degree() > 0
        for s in ws
        for _, r in s
        for v in r.values()
    )
    assert not depends and len(flags) == 8
    assert dims == [[20] * 4, [38] * 4, [72] * 4]
    print("Constant trace limit", job, "verified", flush=True)


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=2) as pool:
        list(pool.map(verify, product((1, -1), (1, -1), range(2))))
        list(pool.map(verify_trace_jump, product((1, -1), (1, -1), (1, -1), range(2))))
