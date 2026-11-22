"""Relative four-cut through degree thirteen on both rank-fifteen surfaces."""

from concurrent.futures import ProcessPoolExecutor
from itertools import product

from flint import fmpz_poly

from verify_rank_fifteen_pencil_trace_disks import local_relations
from verify_rank_fifteen_relative_trace_centers import (
    C,
    POLY_ONE,
    PRIMITIVE,
    ZERO,
    Coefficient,
    check_cut,
)

NEGATED_PRIMITIVE = C * C + C + 1
POLY_ZERO = fmpz_poly([])


class CommonDenominatorJet:
    """A polynomial in the uniformizer with one Q[c] denominator for all terms."""

    __slots__ = ("coeffs", "den")
    precision = 25

    def __init__(self, coeffs, den=POLY_ONE, *, reduce=True):
        if not den:
            raise ZeroDivisionError
        coeffs = coeffs[: self.precision]
        while coeffs and not coeffs[-1]:
            coeffs.pop()
        if not coeffs:
            den = POLY_ONE
        elif reduce and den != POLY_ONE:
            common = den
            for coeff in reversed(coeffs):
                if common == POLY_ONE:
                    break
                if coeff:
                    common = common.gcd(coeff)
            if den.leading_coefficient() < 0:
                common = -common
            if common != POLY_ONE:
                coeffs = [coeff // common for coeff in coeffs]
                den = den // common
        self.coeffs, self.den = coeffs, den

    @classmethod
    def constant(cls, value):
        value = Coefficient.constant(value)
        return cls([value.num], value.den, reduce=False)

    @classmethod
    def from_fractions(cls, coeffs):
        den = POLY_ONE
        for coeff in coeffs:
            den = den * (coeff.den // den.gcd(coeff.den))
        return cls([coeff.num * (den // coeff.den) for coeff in coeffs], den)

    def __bool__(self):
        return bool(self.coeffs)

    def __neg__(self):
        return type(self)([-a for a in self.coeffs], self.den, reduce=False)

    def __add__(self, other):
        cls = type(self)
        if not isinstance(other, cls):
            other = cls.constant(other)
        if not self:
            return other
        if not other:
            return self
        left, right = self.coeffs, other.coeffs
        size = min(self.precision, max(len(left), len(right)))
        if self.den == other.den:
            return cls(
                [
                    (left[i] if i < len(left) else POLY_ZERO)
                    + (right[i] if i < len(right) else POLY_ZERO)
                    for i in range(size)
                ],
                self.den,
            )
        common = self.den.gcd(other.den)
        a, b = other.den // common, self.den // common
        return cls(
            [
                (left[i] * a if i < len(left) else POLY_ZERO)
                + (right[i] * b if i < len(right) else POLY_ZERO)
                for i in range(size)
            ],
            self.den * a,
        )

    __radd__ = __add__

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        cls = type(self)
        if not isinstance(other, cls):
            other = cls.constant(other)
        if not self or not other:
            return cls([])
        if len(self.coeffs) == 1 and self.coeffs[0] == self.den:
            return other
        if len(other.coeffs) == 1 and other.coeffs[0] == other.den:
            return self
        size = min(self.precision, len(self.coeffs) + len(other.coeffs) - 1)
        coeffs = [POLY_ZERO] * size
        for i, a in enumerate(self.coeffs[:size]):
            if a:
                for j, b in enumerate(other.coeffs[: size - i]):
                    if b:
                        coeffs[i + j] = coeffs[i + j] + a * b
        return cls(coeffs, self.den * other.den)

    __rmul__ = __mul__

    def unit(self):
        assert self.den % PRIMITIVE and self.den % NEGATED_PRIMITIVE
        return (
            bool(self)
            and bool(self.coeffs[0] % PRIMITIVE)
            and bool(self.coeffs[0] % NEGATED_PRIMITIVE)
        )

    def inverse(self):
        assert self.unit()
        if len(self.coeffs) == 1:
            return type(self)([self.den], self.coeffs[0])
        coeffs = [Coefficient(a, self.den) for a in self.coeffs[: self.precision]]
        first = coeffs[0].inverse()
        values = [first]
        for n in range(1, self.precision):
            rest = sum(
                (coeffs[k] * values[n - k] for k in range(1, min(len(coeffs), n + 1))),
                ZERO,
            )
            values.append(-first * rest)
        return type(self).from_fractions(values)

    def __truediv__(self, other):
        return self * other.inverse()

    def divide_by_parameter(self):
        assert not self or not self.coeffs[0]
        return type(self)(self.coeffs[1:], self.den, reduce=False)


def certify(job):
    sign_w, color = job
    series = CommonDenominatorJet
    series.precision = 25
    # n=(c²-1)t. This is a change by a unit at all four chosen points.
    relations = local_relations(
        series,
        1,
        1,
        sign_w,
        color,
        level=series.constant(Coefficient(C)),
        parameter=series([POLY_ZERO, C * C - 1]),
    )
    check_cut(relations, series, job, 13)


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=3) as pool:
        list(pool.map(certify, product((1, -1), range(2))))
