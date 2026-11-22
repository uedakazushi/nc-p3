"""Exact arithmetic on Y² = 2(t-1)(t+1)(t+3) over Q(i) and the four relation spaces."""

from functools import cache
from operator import index

import sympy as sp
from flint import fmpz_poly

from verify_fixed_geometric_quartic_locus import a, e, f, symbolic_blocks
from verify_rank_fourteen_descent import frame, slopes
from verify_uniform_fixed_quartic_residuals import b
from verify_uniform_fixed_quartic_support import l, z

ZERO = fmpz_poly([])
ONE = fmpz_poly([1])
VARIABLE = fmpz_poly([0, 1])
CUBIC = 2 * (VARIABLE - 1) * (VARIABLE + 1) * (VARIABLE + 3)


def gaussian_product(a, b, c, d):
    return a * c - b * d, a * d + b * c


class CurveElement:
    """(a + b*i + c*Y + d*i*Y)/q with primitive integral polynomial coordinates."""

    __slots__ = ("coordinates", "denominator")

    def __init__(self, coordinates, denominator=ONE):
        if not denominator:
            raise ZeroDivisionError
        if not any(coordinates):
            self.coordinates = (ZERO,) * 4
            self.denominator = ONE
            return
        if denominator != ONE:
            common = denominator
            for coefficient in coordinates:
                if common == ONE:
                    break
                common = common.gcd(coefficient)
            if denominator.leading_coefficient() < 0:
                common = -common
            if common != ONE:
                coordinates = tuple(c // common for c in coordinates)
                denominator = denominator // common
        self.coordinates, self.denominator = coordinates, denominator

    @classmethod
    def lift(cls, value):
        if isinstance(value, cls):
            return value
        poly = value if isinstance(value, fmpz_poly) else fmpz_poly([index(value)])
        return cls((poly, ZERO, ZERO, ZERO))

    def __bool__(self):
        return any(self.coordinates)

    def __eq__(self, other):
        other = self.lift(other)
        return (
            self.coordinates == other.coordinates
            and self.denominator == other.denominator
        )

    def __neg__(self):
        return CurveElement(tuple(-v for v in self.coordinates), self.denominator)

    def __add__(self, other):
        other = self.lift(other)
        if not other:
            return self
        if not self:
            return other
        if self.denominator == other.denominator:
            return CurveElement(
                tuple(a + b for a, b in zip(self.coordinates, other.coordinates)),
                self.denominator,
            )
        common = self.denominator.gcd(other.denominator)
        left, right = self.denominator // common, other.denominator // common
        return CurveElement(
            tuple(
                a * right + b * left
                for a, b in zip(self.coordinates, other.coordinates)
            ),
            self.denominator * right,
        )

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other) + (-self)

    def __mul__(self, other):
        other = self.lift(other)
        if not self or not other:
            return zero
        a, b, c, d = self.coordinates
        e, f, g, h = other.coordinates
        x, y = gaussian_product(a, b, e, f)
        u, v = gaussian_product(c, d, g, h)
        m, n = gaussian_product(a + c, b + d, e + g, f + h)
        return CurveElement(
            (x + CUBIC * u, y + CUBIC * v, m - x - u, n - y - v),
            self.denominator * other.denominator,
        )

    __rmul__ = __mul__

    def quadratic_norm(self):
        a, b, c, d = self.coordinates
        return a * a - b * b - CUBIC * (c * c - d * d), 2 * (a * b - CUBIC * c * d)

    def norm_numerator(self):
        x, y = self.quadratic_norm()
        return x * x + y * y

    def inverse(self):
        a, b, c, d = self.coordinates
        x, y = self.quadratic_norm()
        m, n = gaussian_product(a, b, x, -y)
        u, v = gaussian_product(-c, -d, x, -y)
        return CurveElement(
            tuple(self.denominator * k for k in (m, n, u, v)), x * x + y * y
        )

    def __truediv__(self, other):
        return self * self.lift(other).inverse()

    def __rtruediv__(self, other):
        return self.lift(other) * self.inverse()

    def __pow__(self, exponent):
        if exponent < 0:
            return self.inverse() ** (-exponent)
        if exponent == 0:
            return one
        half = self ** (exponent // 2)
        square = half * half
        return square * self if exponent % 2 else square


zero = CurveElement.lift(0)
one = CurveElement.lift(1)
parameter = CurveElement.lift(VARIABLE)
imaginary = CurveElement((ZERO, ONE, ZERO, ZERO))
ordinate = CurveElement((ZERO, ZERO, ONE, ZERO))


def relations():
    r, lam = slopes(parameter)
    s = ordinate / ((parameter + 1) * (parameter + 3))
    assert imaginary**2 == -one and s * s == lam * r
    X, bb, ee, ff = frame(s, imaginary)
    values = {z: r * r, l: lam, a: r * X, b: bb, e: r * ee, f: ff}

    @cache
    def evaluate(expression):
        if expression in values:
            return values[expression]
        if expression.is_Rational:
            return CurveElement.lift(int(expression.p)) / int(expression.q)
        if expression.is_Pow:
            return evaluate(expression.base) ** index(expression.exp)
        if expression.is_Add:
            return sum(map(evaluate, expression.args), zero)
        if expression.is_Mul:
            result = one
            for factor in expression.args:
                result *= evaluate(factor)
            return result
        raise ValueError(expression)

    sides = symbolic_blocks(1, z, l, b)
    return [
        [
            (h, {(j, j ^ h): evaluate(v) for j, v in enumerate(row) if v != 0})
            for h, (rows, _) in enumerate(side)
            for row in rows
        ]
        for side in sides
    ]


def irreducible_factors(poly):
    return {tuple(map(int, g)) for g, _ in poly.factor()[1]}


def excluded_factors():
    r, lam = slopes(parameter)
    zz = r * r
    factors = set()
    for value in (
        parameter - 1,
        parameter + 1,
        parameter + 3,
        zz,
        lam,
        zz * zz - 1,
        lam * lam * zz - 1,
        lam * lam - 1,
        lam * lam * zz * zz - 1,
    ):
        factors.update(irreducible_factors(value.coordinates[0]))
        factors.update(irreducible_factors(value.denominator))
    return factors
