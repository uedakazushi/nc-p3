"""Prove the root exchange and separate the two classes off the branch divisor."""

from itertools import product

import sympy as sp
from flint import fmpq_mpoly_ctx

from flint_rational_function import rational_function_field
from verify_rank_fifteen_coordinate_collisions import oriented_cubic
from verify_rank_fifteen_endpoint_orbits import gram_component, multiply


def coefficients(z, lam):
    return (lam + 1) * (lam * z + 1) / (4 * lam * z), -((lam - 1) * (lam * z - 1)) / (
        4 * lam * z
    )


def operator(z, lam, n, unit):
    k, m = coefficients(z, lam)
    component = gram_component((z, lam, 1 + n, k + m / n, 1 - n, k - m / n), unit)
    right = (unit, unit, -unit, -unit)
    left = (unit, unit, -unit, -1 / z)
    return [
        [
            left[b] * component(i, j, b, t) * right[t]
            for b, t in product(range(4), repeat=2)
        ]
        for i, j in product(range(4), repeat=2)
    ]


def prove_exchange():
    names = ("s", "d", "n", "ii")
    field = rational_function_field(names)
    s, d, n, ii = map(field.from_sympy, sp.symbols(names))
    r, lam = s * s, d * d
    z = r * r
    _, m = coefficients(z, lam)
    source = operator(z, lam, n, field.one)
    target = operator(z, 1 / (z * lam), r * m / n, field.one)
    scale = (field.one, field.one, 1 / (d * s), 1 / (d * s))
    a0 = (field.one, field.one, 1 / r, r)
    a1 = (1 / r, r, field.one, field.one)

    def local(row, col, weights, sign):
        if row == col:
            return scale[row]
        if row ^ col == 1:
            return scale[row] * sign * ii * weights[row]
        return field.zero

    basis = [
        [
            local(a, c, a0, -1) * local(b, t, a1, 1)
            for c, t in product(range(4), repeat=2)
        ]
        for a, b in product(range(4), repeat=2)
    ]
    modulus = (ii * ii + 1).n
    counts = [0, 0]
    for a, b in product(range(16), repeat=2):
        diff = sum(
            (
                target[a][j] * basis[j][b] - basis[a][j] * source[j][b]
                for j in range(16)
            ),
            field.zero,
        )
        assert not diff.n % modulus
        counts[bool(diff)] += 1
    assert counts == [224, 32]


def prove_separation():
    names = ("n", "z", "lam")
    field = rational_function_field(names)
    n, z, lam = map(field.from_sympy, sp.symbols(names))
    k, m = coefficients(z, lam)
    middle = 1 - z * k * k
    quartic = (n * n * n * n + middle * n * n + z * m * m).n
    discriminant = (middle * middle - 4 * z * m * m).n
    tensor = operator(z, lam, n, field.one)
    blocks = [
        [[tensor[4 * i + (u ^ h)][4 * (i ^ h) + u] for u in range(4)] for i in range(4)]
        for h in range(4)
    ]
    ctx = fmpq_mpoly_ctx.get(names)
    denominator = ctx.constant(1)
    for block in blocks:
        for row in block:
            for x in row:
                if x:
                    denominator = denominator // denominator.gcd(x.d) * x.d
    integral = [
        [
            [x.n * (denominator // x.d) if x else ctx.constant(0) for x in row]
            for row in block
        ]
        for block in blocks
    ]
    adjoint = [
        [[a[k ^ h][j ^ h] for k in range(4)] for j in range(4)]
        for h, a in enumerate(integral)
    ]
    inputs = [multiply(u, v) for u, v in zip(adjoint, integral)]
    outputs = [multiply(v, u) for u, v in zip(adjoint, integral)]
    difference = oriented_cubic(inputs) - oriented_cubic(outputs)
    n, z, lam = ctx.gens()
    # Clear the invertible leading coefficient before polynomial division.
    remainder = (difference * (16 * z * lam**2) ** difference.degrees()[0]) % quartic
    assert remainder and remainder.degrees()[0] < 4
    for j, x in enumerate((n, z, lam)):
        remainder = remainder // x ** min(powers[j] for powers in remainder.monoms())
    resultant = remainder.resultant(quartic, 0)
    expected = (
        discriminant**2
        * (z - 1) ** 4
        * z**28
        * (z * lam**2 - 1) ** 52
        * lam**56
        * (lam - 1) ** 60
        * (z * lam - 1) ** 60
    )
    assert (
        resultant
        and resultant / resultant.leading_coefficient()
        == expected / expected.leading_coefficient()
    )


if __name__ == "__main__":
    prove_exchange()
    prove_separation()
    print(
        "256 root-exchange identities; the ordinary separation resultant is a unit times delta squared"
    )
