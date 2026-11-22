"""Indexed parameter symmetries and the generic fibers on each coordinate-end curve."""

from itertools import product

import numpy as np
import sympy as sp
from flint import nmod

from verify_rank_fifteen_coordinate_ends import (
    EllipticChart,
    GaussianChart,
    ONE,
    VARIABLE,
    ZERO,
    elliptic_frame,
    gaussian_frame,
)
from verify_rank_fifteen_word_symmetries import prove_tensor_identities


def multiply(A, B):
    columns = list(zip(*B))
    return [[sum(a * b for a, b in zip(row, col)) for col in columns] for row in A]


def trace(M):
    return sum(row[i] for i, row in enumerate(M))


def trace_blocks(frame):
    z, lam, a, b, e, f = frame
    unit = z**0
    X = [[unit, 1 / lam], [1 / lam, z]]
    Y = [[unit, lam], [lam * z, unit]]

    def thick(a, b, scale):
        block = [[a, b], [z * b, a]]
        delta = a * a - z * b * b
        inverse = [[a / delta, -b / delta], [-z * b / delta, a / delta]]
        end = multiply(multiply(Y, inverse), X)
        return [block[j] + [v * scale for v in X[j]] for j in range(2)] + [
            [v / scale for v in Y[j]] + end[j] for j in range(2)
        ]

    blocks = [thick(a, b, unit), list(zip(*thick(e, f, lam)))]

    def component(i, j, k, t):
        h = i ^ j
        if h != k ^ t:
            return 0
        if h < 2:
            return blocks[h][i][k]
        return (z if i == 3 and j == 1 else unit) * (
            z if k == 3 and (t ^ 1) == 1 else unit
        )

    right = (unit, unit, -unit, -unit)
    left = (unit, unit, -unit, -1 / z)
    for i in range(4):
        assert not sum(component(i, j, i, j) * right[j] for j in range(4))
        assert not sum(left[j] * component(j, i, j, i) for j in range(4))
    # T(E[u,u^h]) in the output coordinates (i,i^h), including its scalar kernel.
    T = [
        [
            [left[i ^ h] * component(i, u ^ h, i ^ h, u) * right[u] for u in range(4)]
            for i in range(4)
        ]
        for h in range(4)
    ]
    S = [
        multiply([[M[k ^ h][j ^ h] for k in range(4)] for j in range(4)], M)
        for h, M in enumerate(T)
    ]
    return T, S


def rational_traces(field, frame, t):
    _, S = trace_blocks(frame(field.lift(VARIABLE), field((ZERO, ZERO, ONE, ZERO))))
    denominator = sum(map(trace, S))
    powers = S
    result = []
    for m in (2, 3):
        powers = [multiply(x, y) for x, y in zip(powers, S)]
        v = sum(map(trace, powers)) / denominator**m
        a, _, b, _ = v.coordinates

        def polynomial(g):
            return sum(int(c) * t**i for i, c in enumerate(g))

        if field is EllipticChart:
            assert not b
        result.append(
            sp.cancel(
                (polynomial(a) + sp.I * polynomial(b)) / polynomial(v.denominator),
                extension=sp.I,
            )
        )
    return result


def fiber(t, source, target):
    equations = []
    for a, b in zip(source, target):
        numerator, denominator = sp.fraction(b)
        row = sp.Poly(numerator - a.subs(t, 2) * denominator, t, extension=sp.I)
        # The projective infinity point does not belong to this fiber.
        assert row.degree() == max(sp.degree(numerator, t), sp.degree(denominator, t))
        equations.append(row)
    return equations[0].gcd(equations[1]).monic().as_expr()


def residue(frame, p):
    T, S = trace_blocks(frame)
    den = sum(map(trace, S))
    assert den

    def triangle(h, g, k, x, y, z):
        return sum(S[h][a][x] * S[g][a ^ h][y] * S[k][a ^ h ^ g][z] for a in range(4))

    cubic = 0
    for i, j, k in product(range(4), repeat=3):
        h, g, e = i ^ j, j ^ k, k ^ i
        cubic += triangle(h, g, e, j, k, i) + triangle(g, h, e, k, j, i)

    # T as an element of End(output) tensor End(input), via the trace pairing.
    entries = np.zeros((4, 4, 4, 4), dtype=np.int64)
    for h, i, j in product(range(4), repeat=3):
        entries[i, i ^ h, j ^ h, j] = int(T[h][i][j])
    upper = (4, 0, 3, 1, 2)
    lower = (2, 4, 3, 0, 1)
    labels = "abcdefghij"
    forms = [
        labels[i] + labels[upper[i]] + labels[5 + i] + labels[5 + lower[i]]
        for i in range(5)
    ]
    assert 4**10 * (p - 1) ** 5 < 2**63
    contraction = int(
        np.einsum(",".join(forms) + "->", *([entries] * 5), optimize=True)
    )
    return int(2 * cubic / den**3), int(nmod(contraction % p, p) ** 2 / den**5)


def verify():
    prove_tensor_identities()
    t = sp.Symbol("t")
    elliptic = rational_traces(EllipticChart, elliptic_frame, t)
    gaussian = rational_traces(GaussianChart, gaussian_frame, t)
    linear = (t - 2) * (t - sp.Rational(1, 2))
    assert (
        sp.expand(
            fiber(t, elliptic, elliptic)
            - linear * (t * t - sp.Rational(110, 9) * t + 1)
        )
        == 0
    )
    assert (
        sp.expand(
            fiber(t, gaussian, gaussian)
            - (t - 2) * (t - 3) * (t + sp.Rational(1, 2)) * (t + sp.Rational(1, 3))
        )
        == 0
    )
    assert fiber(t, elliptic, gaussian) == 1
    conjugate = [value.xreplace({sp.I: -sp.I}) for value in gaussian]
    assert (
        sp.expand(
            fiber(t, gaussian, conjugate)
            - (t + 2) * (t + 3) * (t - sp.Rational(1, 2)) * (t - sp.Rational(1, 3))
        )
        == 0
    )

    for x, y, expected in (
        (2, 47, (65, 99)),
        (2, 56, (63, 99)),
        (52, 14, (65, 99)),
        (52, 89, (63, 99)),
        (25, 1, (63, 71)),
        (25, 102, (65, 71)),
        (33, 44, (63, 71)),
        (33, 59, (65, 71)),
    ):
        p = 103
        assert (y * y + x**3 - 14 * x * x + x) % p == 0
        assert residue(elliptic_frame(nmod(x, p), nmod(y, p)), p) == expected
    for x, expected in ((2, (93, 9)), (3, (98, 9)), (67, (98, 9)), (50, (93, 9))):
        assert residue(gaussian_frame(nmod(x, 101), nmod(10, 101)), 101) == expected
    for sign, rows in (
        (1, ((2, 95), (3, 99), (36, 99), (54, 95))),
        (-1, ((107, 95), (106, 99), (73, 99), (55, 95))),
    ):
        for x, cubic in rows:
            values = gaussian_frame(nmod(x, 109), nmod(sign * 33, 109))
            z, lam, a, b, e, f = values
            assert z * (b * b + f * f) - 2 * a * e
            assert residue(values, 109) == (cubic, 46 if sign == 1 else 63)
    print(
        "Both multiplier transports hold. Each normalized endpoint component has generic indexed fiber of size two."
    )


if __name__ == "__main__":
    verify()
