"""Identify the companion matrix with the polar map and account for all sixteen contact planes."""

from itertools import combinations

from flint import fmpq, fmpq_mpoly_ctx
import sympy as sp

from verify_klein_elliptic_completion import coeffs
from verify_kummer_discriminant import discriminant_form


def polarity(z, l):
    u = z * (l * l + 1)
    v = l * (z + 1)
    w = l * l * z * z + 1
    return sp.Matrix(
        [
            [u, v, z * v, w],
            [v, l * l + 1, w, v],
            [z * v, w, z * u, z * v],
            [w, v, z * v, u],
        ]
    )


def check():
    variables, form = discriminant_form()
    a, b, c, d, z, l = variables
    t = polarity(z, l)
    assert sp.factor(t.det() - (z * z - 1) ** 2 * (l * l * z - 1) ** 4) == 0
    co, p, q = coeffs(*variables)
    f, F, g, G, e, E, h, H = co
    companion = [
        (1 + e * (H - E)) / (2 * H),
        (f * G - g * F - 1) / (2 * G),
        (g * F - f * G - 1) / (2 * F),
        (1 + e * (E - H)) / (2 * E),
    ]
    ctx = fmpq_mpoly_ctx.get(tuple(map(str, variables)))

    def polynomial(expr):
        return ctx.from_dict(
            {
                key: fmpq(int(c.p), int(c.q))
                for key, c in sp.Poly(expr, *variables).terms()
            }
        )

    fractions = [
        list(map(polynomial, sp.fraction(sp.cancel(entry)))) for entry in companion
    ]
    quartic = polynomial(form)
    gradient = [quartic.derivative(i) for i in range(4)]
    polar = [sum(polynomial(t[i, j]) * gradient[j] for j in range(4)) for i in range(4)]
    for i, j in combinations(range(4), 2):
        ni, di = fractions[i]
        nj, dj = fractions[j]
        _, remainder = divmod(ni * dj * polar[j] - nj * di * polar[i], quartic)
        assert not remainder
    print("All six projective polar identities hold in the quartic coordinate ring.")

    images = polar + list(ctx.gens()[4:])
    assert not divmod(quartic.compose(*images), quartic)[1]
    swapped = [sp.Rational(36, 19), sp.Rational(36, 37), sp.Rational(36, 5), -12]
    _, pp, qq = coeffs(*swapped, sp.Integer(3), sp.Integer(2))
    assert pp[0] * qq[0] != 0
    assert pp[1] ** 2 == 4 * pp[0] * pp[2]
    assert qq[1] ** 2 == 4 * qq[0] * qq[2]
    print(
        "The polar map preserves the quartic; the swapped rational sample has two genuine double roots."
    )

    planes = [
        b - c,
        b * z * z - c,
        a * z - d,
        -a + d * z,
        b - d * l,
        b * l * z - d,
        -c + d * l * z,
        c * l - d,
        -a + b * l * z,
        a * l - b,
        a - c * l,
        a * l * z - c,
    ]
    planes += [sp.factor(co[i]).as_numer_denom()[0] for i in (1, 3, 5, 7)]
    excluded = polynomial(
        z * l * (z * z - 1) * (l * l * z - 1) * (l * l - 1) * (l * l * z * z - 1)
    )
    variables_ring = list(ctx.gens())

    def unit(value):
        return bool(value) and all(
            not divmod(excluded, f)[1] for f, _ in value.factor()[1]
        )

    nodes = []
    for i, plane in enumerate(planes):
        normal = sp.Matrix([sp.diff(plane, v) for v in variables[:4]])
        normals = list(map(polynomial, normal))
        pivot = next(j for j, h in enumerate(normals) if unit(h))
        free = [j for j in range(4) if j != pivot]
        vector = [normals[pivot] * q for q in variables_ring[:4]]
        vector[pivot] = -sum(normals[j] * variables_ring[j] for j in free)
        restricted = quartic.compose(*(vector + variables_ring[4:]))
        factors = restricted.factor()[1]
        curves = [(f, power) for f, power in factors if any(f.degrees()[:4])]
        assert len(curves) == 1 and curves[0][1] == 2
        assert all(unit(f) for f, power in factors if not any(f.degrees()[:4]))
        conic = curves[0][0]
        h = [[conic.derivative(j).derivative(k) for k in free] for j in free]
        det = sum(
            h[0][j]
            * (
                h[1][(j + 1) % 3] * h[2][(j + 2) % 3]
                - h[1][(j + 2) % 3] * h[2][(j + 1) % 3]
            )
            for j in range(3)
        )
        assert unit(det)
        node = [sp.factor(v) for v in t * normal]
        scale = next(v for v in node if v != 0)
        node = [sp.factor(v / scale) for v in node]
        nodes.append(node)
        determinant = sp.factor(node[0] * node[3] - node[1] * node[2])
        if 4 <= i < 12:
            _, first, second = coeffs(*node, z, l)
            assert determinant != 0 and all(
                sp.cancel(v) == 0 for v in (*first, *second)
            )
        else:
            assert determinant == 0
        # A singular point for every coefficient, not only for the test pair.
        # Substitute rational polynomial coordinates after clearing a common denominator.
        common = sp.lcm_list([sp.denom(v) for v in node])
        substitution = [polynomial(sp.cancel(common * v)) for v in node] + list(
            ctx.gens()[4:]
        )
        assert all(not v.compose(*substitution) for v in gradient)
    assert len({tuple(n.subs({z: 3, l: 2}) for n in node) for node in nodes}) == 16
    print(
        "Eight nodes have invertible A0 and simultaneous-zero polynomials; eight have determinant zero."
    )


if __name__ == "__main__":
    check()
