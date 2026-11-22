"""Arithmetic in Q(i, sqrt(3))(t,Y), with Y²=t⁸+14t⁴+1."""

from functools import cache

from rank_fifteen_octic import octic_frame
from rank_fourteen_function_field import (
    CurveElement,
    VARIABLE,
    ZERO,
    ONE,
    evaluate_relations,
)
from verify_fixed_geometric_quartic_locus import a, e, f
from verify_uniform_fixed_quartic_residuals import b
from verify_uniform_fixed_quartic_support import l, z


class OcticFunction(CurveElement):
    __slots__ = ()
    radicand = VARIABLE**8 + 14 * VARIABLE**4 + 1


class Coefficient:
    """a+b*sqrt(3) with a and b in Q(i)(t,Y)."""

    __slots__ = ("a", "b", "reciprocal")

    def __init__(self, a, b=None):
        self.a = a
        self.b = b if b else None
        self.reciprocal = None

    @classmethod
    def lift(cls, value):
        if isinstance(value, cls):
            return value
        return cls(OcticFunction.lift(value))

    def __bool__(self):
        return bool(self.a) or bool(self.b)

    def __neg__(self):
        return Coefficient(-self.a, -self.b if self.b else None)

    def __add__(self, other):
        other = self.lift(other)
        radical = self.b + other.b if self.b and other.b else self.b or other.b
        return Coefficient(self.a + other.a, radical)

    __radd__ = __add__

    def __sub__(self, other):
        return self + -self.lift(other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        other = self.lift(other)
        if self.b and other.b:
            ac = self.a * other.a
            bd = self.b * other.b
            mixed = (self.a + self.b) * (other.a + other.b) - ac - bd
            return Coefficient(ac + 3 * bd, mixed)
        radical = self.b * other.a if self.b else self.a * other.b if other.b else None
        return Coefficient(self.a * other.a, radical)

    __rmul__ = __mul__

    def norm(self):
        square = self.a * self.a
        return square - 3 * self.b * self.b if self.b else square

    def inverse(self):
        if self.reciprocal is None:
            if self.b:
                scale = self.norm().inverse()
                self.reciprocal = Coefficient(self.a * scale, -self.b * scale)
            else:
                self.reciprocal = Coefficient(self.a.inverse())
        return self.reciprocal

    def __truediv__(self, other):
        return self * self.lift(other).inverse()

    def __rtruediv__(self, other):
        return self.inverse() * other

    def __pow__(self, exponent):
        if exponent < 0:
            return self.inverse() ** (-exponent)
        if exponent == 0:
            return self.lift(1)
        half = self ** (exponent // 2)
        square = half * half
        return square * self if exponent % 2 else square


UNIT = Coefficient.lift(1)


@cache
def windows(color):
    t = Coefficient.lift(OcticFunction.lift(VARIABLE))
    y = Coefficient.lift(OcticFunction((ZERO, ZERO, ONE, ZERO)))
    i = OcticFunction((ZERO, ONE, ZERO, ZERO))
    j = Coefficient(i, OcticFunction.lift(1)) / 2
    assert not (j**4 - j * j + 1) and not (y * y - t**8 - 14 * t**4 - 1)
    frame = octic_frame(t, y, j, color)
    return evaluate_relations(dict(zip((z, l, a, b, e, f), frame)), UNIT)
