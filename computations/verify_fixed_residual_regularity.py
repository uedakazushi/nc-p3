"""Central-quartic certificates for the endpoint, companion and octic windows at (3,2)."""

import sympy as sp

from verify_fixed_branch_regular import check_window, red
from verify_fixed_exchange_isogeny import exchange
from verify_fixed_quintic_locus import a, check
from verify_uniform_fixed_quartic_residuals import companion, r, l

z, slope = sp.Integer(3), sp.Integer(2)


def certify(aa, bb, ee, ff, prime, label):
    A = sp.Matrix([[aa, bb], [z * bb, aa]])
    B = sp.Matrix([[ee, ff], [z * ff, ee]])
    assert all(
        red(v, prime)
        for v in (
            A.det(),
            B.det(),
            A.det() + B.det(),
            aa * aa - z * ff * ff,
            ee * ee - z * bb * bb,
        )
    )
    check_window(A, B, z, slope, label, prime)


def check_families():
    prime = 193
    roots = [v for v in range(prime) if v * v % prime == 3]
    assert roots == [14, 179]
    for kind in ("zero", "infinity", 1, -1):
        values = companion(kind)
        for rr in roots[:1] if isinstance(kind, str) else roots:
            point = [red(v.subs({r: rr, l: slope}), prime) for v in values]
            certify(*point, prime, (kind, rr))


def check_octics():
    prime = 1009
    residual = check(1)
    assert len(residual) == len({P for P, _ in residual}) == 2
    print(len(residual), "distinct octic closed points", flush=True)
    reductions = set()
    for P, U in sorted(residual, key=lambda pair: str(pair[0])):
        poly = sp.Poly(P, a)
        assert poly.degree() == 8 and poly.is_irreducible
        assert all(c.q % prime for c in poly.all_coeffs() + sp.Poly(U, a).all_coeffs())
        aa = next(
            v
            for v in range(prime)
            if red(poly.eval(v), prime) == 0 and red(poly.diff().eval(v), prime)
        )
        uv = red(U.subs(a, aa), prime)
        reductions.add((aa, uv))
        N, D = exchange(1, z, slope) * sp.Matrix([aa, 1])
        certify(aa, 1, N * uv, D * uv, prime, ("octic", aa, uv))
    assert reductions == {(60, 14), (33, 982)}


if __name__ == "__main__":
    check_families()
    check_octics()
