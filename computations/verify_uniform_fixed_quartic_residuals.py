"""Remove the companion sheets and certify both coordinate fibers at every ordinary slope."""

import sympy as sp
from verify_fixed_geometric_quartic_locus import minor_indices, symbolic_blocks, a, e, f
from verify_fixed_exchange_isogeny import exchange
from verify_fixed_line_elliptic_curves import bad_sections
from verify_fixed_quartic_ranks import projected_entries
from verify_klein_elliptic_completion import coeffs
from verify_opposite_fixed_determinants import section
from verify_uniform_fixed_quartic_support import (
    cleared_rows,
    determinant,
    polynomial,
    z,
    l,
    u,
)

r, b = sp.symbols("r b")


def signed_equation(sign):
    N, D = exchange(sign, z, l) * sp.Matrix([a, b])
    co, _, _ = coeffs(a, b, sign * z * b, sign * a, z, l)
    delta = a * a - z * b * b
    common = l * (z + sign) * delta
    first = sp.cancel(
        common * (D * u - co[0] - sign * co[1] * (a * a - z * D**2 * u * u)) / D
    )
    second = sp.cancel(
        common * (N * u - co[4] - sign * co[5] * (N * N * u * u - z * b * b)) / N
    )
    assert sp.expand(first - second) == 0 and sp.denom(first) == 1
    det = exchange(sign, z, l).det()
    assert sp.expand(N * N - z * D * D + det * delta) == 0
    if sign == 1:
        eq = sp.factor(
            first.subs({z: r * r, u: 1 / (r * (l * l * r * r - 1))}, simultaneous=True)
        )
        lines = sp.S.One
        for eps in (1, -1):
            A, _ = section(r, l, eps)
            lines *= A[0, 1] * a - A[0, 0] * b
        ratio = sp.factor(eq / lines)
        assert ratio != 0 and ratio.free_symbols <= {r, l}
        # Ordinary factors only, so the two lines exhaust even special slopes.
        allowed = ordinary_factors()
        assert scalar_factors(ratio) <= allowed


def ordinary_factors():
    return set().union(
        *(
            scalar_factors(w)
            for w in (
                r,
                l,
                r**2 - 1,
                r**2 + 1,
                l * l * r * r - 1,
                l * l - 1,
                l * l * r**4 - 1,
            )
        )
    )


def scalar_factors(value):
    # Reuse the dense determinant's polynomial context, with its z slot now holding r.
    assert value != 0
    return {
        str(g)
        for part in sp.fraction(sp.factor(value.subs(r, z)))
        for g, _ in polynomial(part).factor()[1]
    }


def count_marked_points():
    # Normalize the scale u to t=sqrt(det M)*u; it distinguishes the two sheets.
    for sign in (1, -1):
        fibers = {}
        for sig, aa, ee, ff in bad_sections(r * r, l):
            if sig != sign:
                continue
            N, D = exchange(sign, r * r, l) * sp.Matrix([aa, 1])
            assert sp.cancel(ee * D - ff * N) == 0
            t = sp.factor(ff * r * (l * l * r * r - 1 if sign == 1 else 1) / D)
            fibers.setdefault(t, []).append(aa)
        assert set(fibers) == {r, -r, 1 / r, -1 / r}
        for first_coords in fibers.values():
            if sign == 1:
                assert len(first_coords) == 1
            else:
                x, y = first_coords
                assert sp.factor(x * y + r * r) == 0
                assert sp.factor((x - y) ** 2 - (l * l * r * r + 1) ** 2 / l**2) == 0
    A, _ = section(r, l, 1)
    B, _ = section(r, l, -1)
    assert scalar_factors(A[0, 0] * B[0, 1] - A[0, 1] * B[0, 0]) <= ordinary_factors()


def companion(kind):
    if kind in ("zero", "infinity"):
        A = sp.Matrix([[0, 1], [r * r, 0]]) if kind == "zero" else sp.eye(2)
        co, _, _ = coeffs(*A, r * r, l)
        N, D = exchange(1, r * r, l) * sp.Matrix([A[0, 0], A[0, 1]])
        if kind == "zero":
            ff = sp.factor(-1 / (r * r * co[1]))
            ee = sp.factor(N * ff / D)
        else:
            ee = sp.factor(1 / co[5])
            ff = sp.factor(D * ee / N)
    else:
        A, B = section(r, l, kind)
        co, _, _ = coeffs(*A, r * r, l)
        ff = sp.factor(-1 / (r * r * co[1]) - B[0, 1])
        ee = sp.factor(1 / co[5] - B[0, 0])
    return A[0, 0], A[0, 1], ee, ff


def check_point(kind, blocks):
    aa, bb, ee, ff = companion(kind)
    det_a = aa * aa - r * r * bb * bb
    det_b = ee * ee - r * r * ff * ff
    allowed = ordinary_factors()
    for value in (
        det_a,
        det_b,
        det_a + det_b,
        aa * aa - r * r * ff * ff,
        ee * ee - r * r * bb * bb,
    ):
        allowed |= scalar_factors(value)
    subs = {a: aa, b: bb, e: ee, f: ff, z: r * r}

    def denominator_ok(den):
        return scalar_factors(den.subs(z, r)) <= allowed

    for i in range(4):
        for h, target in enumerate((17, 16, 16, 16)):
            original = sp.Matrix(projected_entries(blocks, i, h))
            cache = {
                v: sp.factor(v.subs(subs, simultaneous=True)).subs(r, z)
                for v in set(original)
            }
            mat = original.applyfunc(cache.get)
            sample = mat.subs({z: 2, l: 3})
            snapshot = [
                [int(v.p) * pow(int(v.q), -1, 307) % 307 for v in row]
                for row in sample.tolist()
            ]
            I, J = minor_indices(mat, snapshot, 1)
            assert len(I) == len(J) == target
            det = determinant(cleared_rows(mat.extract(I, J), denominator_ok))
            assert det and all(str(g) in allowed for g, _ in det.factor()[1])
    print(
        kind,
        "all 16 ranks persist when the chart and determinant sum are nonzero",
        flush=True,
    )


if __name__ == "__main__":
    count_marked_points()
    for sign in (1, -1):
        signed_equation(sign)
    blocks = symbolic_blocks(1, z, l, b)
    for kind in ("zero", "infinity", 1, -1):
        check_point(kind, blocks)
