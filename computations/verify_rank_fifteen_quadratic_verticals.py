"""All quadratic horizontal residues at arbitrary geometric levels in characteristic three."""

from itertools import product
from math import gcd, lcm

import sympy as sp
from flint import fq_default_ctx, fq_default_poly_ctx, nmod_poly

from verify_rank_fifteen_genus_two_unit_line import (
    GenusTwo,
    Recorded,
    ONE,
    ZERO,
    T,
    central_divisor,
    full_fibers,
    unit_model,
)
from verify_rank_fifteen_unit_quartic_towers import encode
from verify_rank_fifteen_unit_residues import verify_unit_point


def rational_surface():
    u, q = sp.symbols("u q")
    E = u**6 - u**2
    S = q - q**3
    Tnum = q * E + u * u * S
    # s=S/E and t=Tnum/E: clear E^3 in t^3-t=u^6(s^3-s).
    identity = Tnum**3 - Tnum * E**2 - u**6 * (S**3 - S * E**2)
    assert sp.Poly(identity, u, q, modulus=3).is_zero
    v, b = sp.symbols("v b")

    def AS(x):
        return x**3 - x

    assert sp.Poly(
        AS(v - b) * AS(v + b) - AS(v) ** 2 + AS(b) ** 2, v, b, modulus=3
    ).is_zero


def field_model(s):
    field = fq_default_ctx(3, 4, "a")
    ring = fq_default_poly_ctx(field)
    one, zero, t = ring.one(), ring.zero(), ring.gen()
    i = (-field.one()).sqrt()
    x = s + i
    c = x.sqrt()
    d = 1 / c**3 - 1 / c
    D = t**3 - t
    F = D**2 - ring([d * d])

    class Curve:
        __slots__ = ("coordinates", "denominator")
        divisor = one

        def __init__(self, a, b=zero, den=one):
            if not den:
                raise ZeroDivisionError
            if not a and not b:
                a, b, den = zero, zero, one
            if den != one:
                g = den.gcd(a).gcd(b)
                if g.degree() > 0:
                    a, b, den = a // g, b // g, den // g
                lc = den.leading_coefficient()
                if lc != 1:
                    scale = ring([1 / lc])
                    a, b, den = a * scale, b * scale, den * scale
            self.coordinates = (a, b)
            self.denominator = den

        @classmethod
        def constant(cls, x):
            if isinstance(x, cls):
                return x
            return cls(x if type(x) == type(t) else ring([x]))

        def __bool__(self):
            return any(self.coordinates)

        def __eq__(self, other):
            other = self.constant(other)
            return (
                self.coordinates == other.coordinates
                and self.denominator == other.denominator
            )

        def __neg__(self):
            return Curve(*(-x for x in self.coordinates), self.denominator)

        def __add__(self, other):
            other = self.constant(other)
            if not other:
                return self
            if not self:
                return other
            a, b = self.coordinates
            u, v = other.coordinates
            g = self.denominator.gcd(other.denominator)
            den1 = self.denominator // g
            den2 = other.denominator // g
            return Curve(
                a * den2 + u * den1, b * den2 + v * den1, self.denominator * den2
            )

        __radd__ = __add__

        def __sub__(self, other):
            return self + (-self.constant(other))

        def __rsub__(self, other):
            return self.constant(other) + (-self)

        def __mul__(self, other):
            other = self.constant(other)
            a, b = self.coordinates
            u, v = other.coordinates
            return Curve(
                a * u + b * v * F, a * v + b * u, self.denominator * other.denominator
            )

        __rmul__ = __mul__

        def inverse(self):
            a, b = self.coordinates
            norm = a * a - b * b * F
            if not norm:
                raise ZeroDivisionError
            f = norm
            while f.degree() > 0:
                g = f.gcd(Curve.divisor)
                if g.degree() == 0:
                    for z, e in f.factor_squarefree()[1]:
                        Curve.divisor *= z
                    break
                f = f // g
            return Curve(a * self.denominator, -b * self.denominator, norm)

        def __truediv__(self, other):
            return self * self.constant(other).inverse()

        def __rtruediv__(self, other):
            return self.inverse() * other

        def __pow__(self, e):
            if e < 0:
                return self.inverse() ** (-e)
            if e == 0:
                return self.constant(1)
            square = self ** (e // 2)
            return square * square * self if e % 2 else square * square

    return (
        Curve,
        ring,
        unit_model(
            Curve.constant(t), Curve.constant(i), Curve.constant(c), Curve(zero, one)
        ),
    )


def descend(ring, div):
    t = ring.gen()
    B = -((t**3 - t) ** 2)
    p = div // div.gcd(B**9 - B)
    reflected = ring([(-1 if j % 2 else 1) * a for j, a in enumerate(p)])
    core = p.gcd(reflected)
    for sign, b in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        # Reduction before composition is valid modulo the transported divisor.
        transported = core.compose(sign * (t - b)).monic()
        core = core.gcd((p % transported).compose_mod(sign * t + b, core))
    assert core.is_squarefree()
    norm = core
    conj = core
    for degree in range(3):
        conj = ring([c**3 for c in conj])
        norm *= conj // conj.gcd(norm)
    coeff = []
    for x in norm:
        parts = x.to_list()
        assert not any(int(c) for c in parts[1:])
        coeff.append(int(parts[0]))
    f = nmod_poly(coeff, 3)
    stable, levels = full_fibers(f)
    assert stable == f.degree()
    return core.degree(), stable, levels


def split_unit_coordinates(t, x, y, c):
    for v, _ in ((t**3 - t) ** 2 - y / x**3).roots():
        k, h = 1 / v, c * (1 / v - v)
        sq = [k * k - 1 + s * 2 * h for s in (1, -1)] + [
            2 * (k + delta + s * h) for delta, s in product((1, -1), repeat=2)
        ]
        if not all(r.is_square() for r in sq):
            continue
        a, b, p, u, q, r = [z.sqrt() for z in sq]
        coords = tuple(map(encode, (k, h, (a + b) / 2, p, u, q, r)))
        return coords
    return None


def certify(g, s=0):
    d = len(g) - 1
    if d <= 2:
        return
    degree = lcm(d, 4 if s else 2)
    for mult in (1, 2, 3, 4, 6, 12):
        field = fq_default_ctx(3, degree * mult, "a")
        ring = fq_default_poly_ctx(field)
        t = ring.gen()
        x = s + (-field.one()).sqrt()
        c = x.sqrt()
        roots = sorted([r for r, e in ring(list(g)).roots()], key=encode)
        seen = set()
        points = []
        for y0 in roots:
            if encode(y0) in seen:
                continue
            r = y0
            while encode(r) not in seen:
                seen.add(encode(r))
                r = r**9
            y = -(x**3) * y0
            coords = split_unit_coordinates(t, x, y, c)
            if coords is not None:
                points.append((y, coords))
        if len(points) != gcd(d, 2):
            continue
        for y, coords in points:
            verify_unit_point(field, x, y, coords, ((0, 1, 0), (0, 0, 1)))
        print(
            "Certified",
            g,
            "in degree",
            degree * mult,
            "on",
            len(points),
            "F9-orbits",
            flush=True,
        )
        return
    raise AssertionError(g)


def verify_line(s):
    if s == 0:
        Recorded.divisor = ONE
        Recorded.record.cache_clear()
        v = GenusTwo.constant(T)
        i = GenusTwo((ZERO, ONE, ZERO, ZERO))
        root = GenusTwo((ZERO, ZERO, ONE, ZERO))
        windows = unit_model(v, i, 1 - i, root)
        divisor = central_divisor(
            windows, GenusTwo.constant(1), Recorded, (1205, 1997, 1898, 2147)
        )
        assert divisor.degree() == 6604
        stable, levels = full_fibers(divisor)
        assert stable == 126 and sorted(len(g) - 1 for g in levels) == [
            2,
            2,
            3,
            4,
            4,
            6,
        ]
    else:
        curve, ring, windows = field_model(s)
        degrees = (3240, 3964, 3281, 2317) if s == 1 else (3321, 3604, 3290, 2843)
        divisor = central_divisor(windows, curve.constant(1), curve, degrees)
        expected = (12491, 54, 72, 3) if s == 1 else (12756, 48, 48, 2)
        stable, prime_stable, levels = descend(ring, divisor)
        assert (divisor.degree(), stable, prime_stable, len(levels)) == expected
        assert all(len(g) == 5 for g in levels)
    for g in levels:
        certify(g, s)
    print(
        "Certified all geometric nonzero levels for residue",
        s,
        "+ i and its conjugate.",
        flush=True,
    )


if __name__ == "__main__":
    rational_surface()
    for s in (0, 1, -1):
        verify_line(s)
