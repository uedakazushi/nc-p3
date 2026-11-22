"""Recover the lifted Pauli parameters and enumerate their projective orbit."""

from collections import Counter
from itertools import product

import sympy as sp
from flint import fq_default_ctx

from flint_rational_function import rational_function_field
from verify_rank_fifteen_pauli_comparison import label_permutations
from verify_rank_fifteen_root_exchange import coefficients, operator


def pattern(a, b, L, M, A, B, C, D):
    return (0, b, A, A, L, M, -C, D, a, a + b, B, B, L, M, -D, C)


def coordinates(v):
    return v[8], v[1], v[4], v[5], v[2], v[10], v[15], v[7]


def quadrics(a, b, L, M, A, B, C, D):
    return (
        L * L + M * M - b * (a + b),
        A * A + B * B - a * (a + b),
        C * C + D * D - a * b,
        L * M + A * B,
        L * M + C * D,
    )


def amplitudes(s, d, n, one):
    t = s * d
    r, w = s * s, t * t
    K = (r + 1 / r + w + 1 / w) / 4
    h = (r + 1 / r - w - 1 / w) / 4
    p, q = s + 1 / s, s - 1 / s
    u, v = t + 1 / t, t - 1 / t
    return pattern(
        one,
        (K - 1) / 2,
        n / 2,
        h / (2 * n),
        (u - p) / 4,
        (u + p) / 4,
        (v - q) / 4,
        (v + q) / 4,
    )


def recover(row):
    a, b, L, M, A, B, C, D = coordinates(row)
    return (B - A + D - C) / a, (B + A + D + C) / a, 2 * L / a


def verify_projections():
    names = ("n", "s", "d")
    field = rational_function_field(names)
    n, s, d = map(field.from_sympy, sp.symbols(names))
    r, lam = s * s, d * d
    z = r * r
    k, m = coefficients(z, lam)
    quartic = (n * n * n * n + (1 - z * k * k) * n * n + z * m * m).n
    tensor = operator(z, lam, n, field.one)
    A = ((field.one, field.one, 1 / r, r), (1 / r, r, field.one, field.one))
    B = ((field.one, 1 / r, d, r * d), (field.one, r, d, r * d))
    expected = amplitudes(s, d, n, field.one)
    for x in range(16):
        shift, phase = x % 4, x // 4
        weights = []
        for side in (0, 1):
            first = (shift & 1) ^ ((shift >> 1) if side else 0)
            second = shift >> 1
            char = (0, 3, 2, 1)[phase] if side else phase
            weights.append(
                [
                    (A[side][i] if first else field.one)
                    * (B[side][i ^ first] if second else field.one)
                    * (-1) ** (((i ^ shift) & char).bit_count())
                    for i in range(4)
                ]
            )
        left, right = weights
        square = left[0] * left[shift] * right[0] * right[shift]
        projection = sum(
            (
                tensor[4 * (i ^ shift) + (j ^ shift)][4 * i + j] * left[i] * right[j]
                for i, j in product(range(4), repeat=2)
            ),
            field.zero,
        ) / (16 * square)
        # Product of the two local square-root normalizations for the B generators.
        if shift & 2:
            projection = projection * s * d
        assert not (projection - expected[x]).n % quartic
    for value in quadrics(*coordinates(expected)):
        assert not value.n % quartic
    assert recover(expected) == (s, s * d, n)


def dot(a, b):
    return (a & b).bit_count() % 2


def product_sign(x, y):
    phase = x // 4
    return dot(phase ^ (0, 3, 2, 1)[phase], y % 4)


def signed_actions():
    for perm in label_permutations():
        exponent = [0]
        for bit in (1, 2, 4, 8):
            root = (product_sign(bit, bit) - product_sign(perm[bit], perm[bit])) % 2
            exponent += [
                (
                    exponent[x]
                    + root
                    + 2 * (product_sign(perm[x], perm[bit]) - product_sign(x, bit))
                )
                % 4
                for x in range(bit)
            ]
        assert all(
            (
                exponent[x]
                + exponent[y]
                - exponent[x ^ y]
                + 2 * (product_sign(perm[x], perm[y]) - product_sign(x, y))
            )
            % 4
            == 0
            for x, y in product(range(16), repeat=2)
        )
        for character in range(16):
            yield (
                perm,
                [(e + 2 * dot(character, x)) % 4 for x, e in enumerate(exponent)],
            )


def template_group():
    variables = sp.symbols("a b L M A B C D")
    source = pattern(*variables)
    ideal = sp.groebner(
        quadrics(*variables), *variables, domain=sp.QQ_I, order="grevlex"
    )
    result = []
    for perm, exp in signed_actions():
        target = tuple(sp.I**e * source[j] for j, e in zip(perm, exp))
        coords = coordinates(target)
        if any(sp.expand(u - v) != 0 for u, v in zip(target, pattern(*coords))):
            continue
        assert all(ideal.reduce(q)[1] == 0 for q in quadrics(*coords))
        result.append((perm, exp))
    assert len(result) == 192 and set(
        Counter(tuple(f) for f, w in result).values()
    ) == {4}
    assert len({tuple(f) for f, w in result}) == 48
    return result


def verify_fiber(group):
    field = fq_default_ctx(101, 2, "v")
    s = field(2).sqrt()
    src = amplitudes(s, field(3), field(46), field.one())
    assert all(q.is_zero() for q in quadrics(*coordinates(src)))
    lifted = set()
    parameters = set()
    multipliers = set()
    for perm, exp in group:
        target = [field(10) ** e * src[j] for j, e in zip(perm, exp)]
        s1, t1, n1 = recover(target)
        z, lam = s1**4, (t1 / s1) ** 2
        k, m = coefficients(z, lam)
        assert not (
            (z * z - 1)
            * (lam * lam - 1)
            * (z * lam * lam - 1)
            * (z * z * lam * lam - 1)
            * (1 - z * k * k)
        ).is_zero()
        assert not ((1 - z * k * k) ** 2 - 4 * z * m * m).is_zero()
        assert all(
            u == v * target[8]
            for u, v in zip(target, amplitudes(s1, t1 / s1, n1, field.one()))
        )
        lifted.add(tuple(map(str, (s1, t1, n1))))
        parameters.add(tuple(map(str, (z, lam, n1))))
        multipliers.add(tuple(map(str, (z, lam))))
    assert (len(lifted), len(parameters), len(multipliers)) == (192, 24, 12)


if __name__ == "__main__":
    verify_projections()
    verify_fiber(template_group())
    print(
        "15 coefficients, five quadrics, 192 template actions, 24 generic first-display parameters and 12 multipliers"
    )
