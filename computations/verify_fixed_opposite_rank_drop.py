"""Determine opposite-rank drops on two rational fibers and certify their quartic orbits."""

import sympy as sp

from verify_fixed_branch_regular import check_window, red
from verify_fixed_exchange_isogeny import exchange
from verify_fixed_geometric_quartic_locus import a
from verify_klein_elliptic_completion import coeffs, gram_twist
from verify_uniform_fixed_quartic_residuals import companion, r, l as slope
from verify_uniform_fixed_quartic_support import u

quartics = [
    932 * a**4 - 8229 * a**3 + 18228 * a**2 - 5265 * a - 576,
    64 * a**4 + 1755 * a**3 - 18228 * a**2 + 24687 * a - 8388,
]


collision = 769 * a**4 + 3744 * a**3 - 10296 * a**2 + 134784 * a + 996624


def factors(value):
    return {g for g, _ in sp.factor_list(value)[1]}


def check_endpoints(z, l):
    for kind in ("zero", "infinity"):
        av, bv, e, f = [
            sp.cancel(v.subs({r: sp.sqrt(z), slope: l})) for v in companion(kind)
        ]
        A = sp.Matrix([[av, bv], [z * bv, av]])
        B = sp.Matrix([[e, f], [z * f, e]])
        assert all(
            sp.simplify(x) != 0
            for x in (A.det(), B.det(), av * av - z * f * f, e * e - z * bv * bv)
        )
        assert gram_twist(A, B, z, l)[1].rank() == 16


def check(sign, z, l, candidates, rank):
    if sign == 1:
        check_endpoints(z, l)
    N, D = exchange(sign, z, l) * sp.Matrix([a, 1])
    A = sp.Matrix([[a, 1], [sign * z, sign * a]])
    B = u * sp.Matrix([[N, D], [sign * z * D, sign * N]])
    co, _, _ = coeffs(*A, z, l)
    eq = sp.factor(
        sign * z * co[1] * (D * u) ** 2 + D * u - co[0] - sign * co[1] * a * a
    )
    (G,) = [g for g in factors(sp.fraction(eq)[0]) if g.has(u)]
    _, opp, _ = gram_twist(A, B, z, l)
    boundary = factors(a * (a * a - z) * N * D)
    found = {}
    for h in range(4):
        indices = [4 * i + (i ^ h) for i in range(4)]
        assert all(
            opp[i, j] == 0 for i in indices for j in set(range(16)) - set(indices)
        )
        block = opp.extract(indices, indices).applyfunc(sp.factor)
        num, den = sp.fraction(sp.factor(block.det(method="domain-ge")))
        assert num and factors(den) <= factors(u * (a * a - z))
        projected = set()
        for g in factors(num):
            norm = sp.resultant(G, g, u) if g.has(u) else g
            assert norm
            projected.update(factors(norm) - boundary)
        expected = {candidates[h]} if sign == 1 and h < 2 else set()
        assert projected == expected
        for P in projected:
            assert sp.Poly(P, a).is_irreducible
            basis = sp.groebner([G, P, num], u, a, domain=sp.QQ).polys
            assert (
                len(basis) == 2
                and basis[1].as_expr() == sp.Poly(P, a).monic().as_expr()
            )
            U = sp.expand(u - basis[0].as_expr())
            assert not U.has(u)
            if P in found:
                assert found[P] == U
                continue
            found[P] = U
            prime = 101
            assert all(c.q % prime for c in sp.Poly(U, a).all_coeffs())
            root = next(
                i
                for i in range(prime)
                if P.subs(a, i) % prime == 0 and sp.diff(P, a).subs(a, i) % prime
            )
            uv = red(U.subs(a, root), prime)
            aa, bb = A.subs(a, root), B.subs({a: root, u: uv})
            assert all(
                red(v, prime)
                for v in (
                    aa.det(),
                    bb.det(),
                    aa[0, 0] ** 2 - z * bb[0, 1] ** 2,
                    bb[0, 0] ** 2 - z,
                )
            )
            check_window(
                aa, bb, z, l, ("quartic orbit", h, root, uv), prime, opposite_rank=rank
            )
    assert set(found) == (set(candidates) if sign == 1 else set())
    print(sign, len(found), "quartic rank-drop orbits", flush=True)


if __name__ == "__main__":
    for sign in (1, -1):
        check(sign, sp.Integer(3), sp.Integer(2), quartics, 15)
    check(1, sp.Integer(36), sp.Rational(-1, 4), [collision, collision], 14)
