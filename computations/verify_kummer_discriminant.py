"""The six tangent branch lines of the discriminant surface and a regular ramified window."""

from functools import cache
from itertools import combinations

import sympy as sp

from verify_klein_elliptic_completion import coeffs, four_relations, gram_twist
from verify_klein_full_rank_configurations import partner
from verify_klein_generic_isomorphisms import framed, number_of_maps
from verify_klein_two_doublets import modular_points
from verify_klein_period_four_center import artinian_quotient, commuting_space, paths


@cache
def discriminant_form():
    a, b, c, d, z, l = sp.symbols("a b c d z l")
    _, p, _ = coeffs(a, b, c, d, z, l)
    denominator = l * (z * z - 1) * (a * d - b * c)
    form = sp.cancel(denominator**2 * (p[1] ** 2 - 4 * p[0] * p[2]))
    return (a, b, c, d, z, l), form


def projection():
    (a, b, c, d, z, l), form = discriminant_form()
    in_a = sp.Poly(form, a)
    assert sp.Poly(form, a, b, c, d).is_homogeneous and in_a.degree() == 2
    lines = (b - c, b - d * l, -c + d * l * z, b * z * z - c, c * l - d, b * l * z - d)
    scalar = 16 * (l * l - 1) * (z - 1) ** 2 * (l * l * z * z - 1)
    assert sp.expand(in_a.discriminant() - scalar * sp.prod(lines)) == 0
    e = z * l + 1 / l
    dual = sp.Matrix(
        [[2, 1 + z * z, e], [1 + z * z, 2 * z * z, z * e], [e, z * e, 2 * z]]
    )
    assert (
        sp.factor(
            dual.det() - 2 * z * (l * l - 1) * (z - 1) ** 2 * (l * l * z * z - 1) / l**2
        )
        == 0
    )
    coeff = sp.hessian(in_a.nth(2), (b, c, d))
    assert (coeff + 2 * l * l * dual.adjugate()).applyfunc(sp.cancel) == sp.zeros(3)
    points = [sp.Matrix([sp.diff(f, t) for t in (b, c, d)]) for f in lines]
    assert all(sp.cancel((v.T * dual * v)[0]) == 0 for v in points)
    # The excluded factors account for every possible pairwise collision.
    excluded = z * l * (z * z - 1) * (l * l * z - 1) * (l * l - 1) * (l * l * z * z - 1)
    for u, v in combinations(points, 2):
        common = sp.factor(sp.gcd_list(list(u.cross(v))))
        assert common != 0 and all(
            sp.rem(excluded, f, z, l) == 0 for f, _ in sp.factor_list(common)[1]
        )
    print(
        "The projection discriminant is six distinct tangents to the computed nonsingular conic."
    )


def regular_sample():
    a, b, c, d, z, l = map(sp.Rational, ("36/17", 0, 9, 4, 3, 2))
    co, p, q = coeffs(a, b, c, d, z, l)
    assert (
        p[0] * q[0] and p[1] ** 2 - 4 * p[0] * p[2] == q[1] ** 2 - 4 * q[0] * q[2] == 0
    )
    x, y = -p[1] / (2 * p[0]), -q[1] / (2 * q[0])
    f0, f1, g0, g1, e0, e1, h0, h1 = co
    first = sp.Matrix([[a, b], [c, d]])
    second = sp.Matrix([[e0 + e1 * y, f0 + f1 * x], [g0 + g1 * x, h0 + h1 * y]])
    assert second == sp.Matrix(
        [[sp.Rational(36, 19), sp.Rational(36, 37)], [sp.Rational(36, 5), -12]]
    )
    assert x * y * first.det() * second.det()
    w, opp, rev = gram_twist(first, second, z, l)
    assert w.rank() == rev.rank() == 6 and opp.rank() == 16
    rules = four_relations(w, rev)
    basis, image = paths(rules, 5)
    for i in range(4):
        assert [sum(len(basis[i, n, h]) for h in range(4)) for n in range(6)] == [
            1,
            4,
            10,
            20,
            35,
            56,
        ]
    assert [commuting_space(basis, image, h, 0).cols for h in range(4)] == [3, 2, 0, 0]
    artinian_quotient(rules, basis, image)
    prime = 34381
    first = modular_points(w, prime)
    second = partner(first, w, False, prime)
    third = partner(second, rev, True, prime)
    fourth = partner(third, w.T, False, prime)
    assert partner(fourth, rev.T, True, prime) == first
    flip = [4 * j + i for i in range(4) for j in range(4)]
    for pts, cut in [(second, rev[flip, :]), (third, w.T), (fourth, rev.T[flip, :])]:
        assert set(modular_points(cut, prime)) == set(pts)
    frames = [framed(pts, prime) for pts in (first, second, third, fourth)]
    assert [[number_of_maps(row, col, prime) for col in frames] for row in frames] == (
        4 * sp.eye(4)
    ).tolist()
    print(
        "Four pairwise inequivalent reduced configurations; all four projective stabilizers have order four."
    )
    print(
        "Genuine double roots, opposite rank 16, and four central parameters at the rational ramification point."
    )


if __name__ == "__main__":
    projection()
    regular_sample()
