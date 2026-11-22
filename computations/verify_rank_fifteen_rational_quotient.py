"""Two rational coordinates for the generic indexed quotient of the rank-fifteen surfaces."""

import sympy as sp

from verify_rank_fifteen_parameter_orbits import (
    coordinates,
    pattern,
    quadrics,
    template_group,
)


def relative_forms(a, b, L, M, A, B, C, D):
    c = a + b
    return (
        2 * (a * a + b * b + c * c),
        4 * a * b * c,
        (b + c) * ((b + c) ** 2 - 9 * a * a),
        4 * L * M,
        16 * (L * L - M * M) * (B * B - A * A) * (D * D - C * C),
    )


def relative_sign(value, reference):
    if sp.expand(value - reference) == 0:
        return 1
    assert sp.expand(value + reference) == 0
    return -1


def verify_group():
    variables = sp.symbols("a b L M A B C D")
    coeffs = pattern(*variables)
    reference = relative_forms(*variables)
    relations = sp.groebner(
        quadrics(*variables), *variables, domain=sp.QQ_I, order="grevlex"
    )
    for perm, exp in template_group():
        row = [sp.I**e * coeffs[j] for j, e in zip(perm, exp)]
        p2, p3, d3, h, theta = relative_forms(*coordinates(row))
        assert (
            sp.expand(p2 - reference[0]) == 0 and sp.expand(theta - reference[4]) == 0
        )
        p_sign = relative_sign(p3, reference[1])
        d_sign = relative_sign(d3, reference[2])
        assert relations.reduce(h - p_sign * d_sign * reference[3])[1] == 0


def verify_equations():
    U, V = sp.symbols("U V")
    I = (U + 3) ** 3 / (U - 1) ** 2
    J = V / (U + 3) ** 2
    G2 = U * V * (U - 9) ** 2 / ((U - 1) ** 2 * (U + 3) ** 2)
    Y2 = ((U - 1) ** 2 - 4 * V) * ((U + V - 1) ** 2 - 4 * U * V) / (U + 3) ** 6
    assert sp.cancel(G2 - J * (I - 27)) == 0
    assert sp.cancel(Y2 - (1 / I - J) ** 2 + 4 * J**3) == 0
    x = V * (U + 3) / (U - 1) ** 2
    y = V**3 / (U - 1) ** 4
    assert sp.cancel(x * x * G2 - x**3 + 27 * y) == 0
    assert sp.cancel(I * I * Y2 - (1 - x) ** 2 + 4 * y) == 0
    first = (U**2 + 14 * U + 16 * V - 7) ** 2 - 16 * U * (3 * U - 1) ** 2
    second = (U**2 - 2 * U + 16 * V + 9) ** 2 - 16 * U * (U - 3) ** 2
    third = (V - 4 * (U + 1)) ** 2 - 64 * U
    elliptic = x**4 - x**3 + 32 * x * x * y - 160 * x * y + 256 * y * y + 125 * y
    rational = x**4 + x**3 + 32 * x * x * y + 256 * y * y - 27 * y
    assert sp.cancel(elliptic + V**3 * (U - V - 2) * first / (U - 1) ** 8) == 0
    assert sp.cancel(rational - V**3 * (U + V) * second / (U - 1) ** 8) == 0
    assert sp.cancel(y - (x - 4) ** 2 + ((U - 1) ** 2 - V) * third / (U - 1) ** 4) == 0
    x, T = sp.symbols("x T")
    c = 2 / (3 * sp.sqrt(3))
    f = (1 - x) ** 2 - 4 * x**3 / 27
    assert sp.expand(27 * f + (x - 3) ** 2 * (4 * x - 3)) == 0
    q = (T + f / T) / 2
    p = (T - f / T) / (2 * c)
    assert sp.cancel(q * q - (1 - x) ** 2 + 4 * (x**3 - p**2) / 27) == 0
    assert sp.expand(q + c * p - T) == 0


def verify_radicals_and_reflection():
    t = sp.symbols("t")
    K = (sp.Rational(17, 4) + t * t + t**-2) / 4
    H = (sp.Rational(17, 4) - t * t - t**-2) / 4
    polys = []
    for e in (1, -1):
        f, g = sp.fraction(sp.factor(K * K - 1 + e * 2 * H))
        assert sp.degree(f, t) == 8 and sp.gcd(f, f.diff(t)) == 1 and sp.gcd(f, g) == 1
        polys.append(f)
    assert sp.gcd(*polys) == 1
    args = sp.symbols("a b L M A B C D")
    a, b, L, M, A, B, C, D = args
    before = pattern(*args)
    after = pattern(a, b, L, -M, -A, B, -C, D)
    for x in range(16):
        shift, phase = x % 4, x // 4
        exponent = (shift & phase).bit_count() + (shift & 2).bit_count()
        assert after[x] == (-1) ** exponent * before[x]


if __name__ == "__main__":
    verify_group()
    verify_equations()
    verify_radicals_and_reflection()
    print(
        "192 exact invariant actions, two quadratic extensions, rational inverse in x,T and partial transpose"
    )
