"""Saturate the trace ±1 germs with variable level at both primitive fibers."""

from concurrent.futures import ProcessPoolExecutor
from itertools import product

from flint import fmpq, fmpz_poly

from verify_rank_fifteen_pencil_trace_disks import (
    local_relations,
    local_central_sections,
)
from verify_rank_fourteen_all_ordinary import DIMENSIONS, saturated_projection
from verify_rank_fourteen_central_cut import quotient_layers


POLY_ONE = fmpz_poly([1])
C = fmpz_poly([0, 1])
PRIMITIVE = C * C - C + 1


class Coefficient:
    """An exact rational function in c, with an invertibility test at Phi_6."""

    __slots__ = ("num", "den")

    def __init__(self, num, den=POLY_ONE):
        if not den:
            raise ZeroDivisionError
        if not num:
            den = POLY_ONE
        elif den != POLY_ONE:
            common = num.gcd(den)
            if den.leading_coefficient() < 0:
                common = -common
            num, den = num // common, den // common
        self.num, self.den = num, den

    @classmethod
    def constant(cls, value):
        if isinstance(value, cls):
            return value
        value = fmpq(value)
        return cls(
            fmpz_poly([int(value.numerator)]), fmpz_poly([int(value.denominator)])
        )

    def __bool__(self):
        return bool(self.num)

    def __neg__(self):
        return Coefficient(-self.num, self.den)

    def __add__(self, other):
        if not self:
            return other
        if not other:
            return self
        if self.den == other.den:
            return Coefficient(self.num + other.num, self.den)
        common = self.den.gcd(other.den)
        left, right = self.den // common, other.den // common
        return Coefficient(self.num * right + other.num * left, self.den * right)

    def __mul__(self, other):
        if not self or not other:
            return ZERO
        return Coefficient(self.num * other.num, self.den * other.den)

    def inverse(self):
        return Coefficient(self.den, self.num)

    def unit(self):
        assert self.den % PRIMITIVE
        return bool(self.num % PRIMITIVE)


ZERO = Coefficient.constant(0)
ONE = Coefficient.constant(1)


class Jet:
    """Truncated series in n over the local ring Q[c]_(Phi_6)."""

    __slots__ = ("coeffs",)
    precision = 0

    def __init__(self, coeffs):
        values = list(coeffs[: self.precision])
        while values and not values[-1]:
            values.pop()
        self.coeffs = values

    @classmethod
    def constant(cls, value):
        return cls([Coefficient.constant(value)])

    def __bool__(self):
        return bool(self.coeffs)

    def __neg__(self):
        return Jet([-v for v in self.coeffs])

    def __add__(self, other):
        if not isinstance(other, Jet):
            other = Jet.constant(other)
        left, right = self.coeffs, other.coeffs
        end = min(len(left), len(right), self.precision)
        tail = left[end:] if len(left) > len(right) else right[end:]
        return Jet([left[i] + right[i] for i in range(end)] + tail)

    __radd__ = __add__

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if not isinstance(other, Jet):
            other = Jet.constant(other)
        if not self or not other:
            return Jet([])
        size = min(self.precision, len(self.coeffs) + len(other.coeffs) - 1)
        values = [ZERO] * size
        for i, a in enumerate(self.coeffs[:size]):
            if a:
                for j, b in enumerate(other.coeffs[: size - i]):
                    if b:
                        values[i + j] = values[i + j] + a * b
        return Jet(values)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self * other.inverse()

    def unit(self):
        return bool(self) and self.coeffs[0].unit()

    def inverse(self):
        assert self.unit()
        first = self.coeffs[0].inverse()
        values = [first]
        if len(self.coeffs) > 1:
            for n in range(1, self.precision):
                rest = sum(
                    (
                        self.coeffs[j] * values[n - j]
                        for j in range(1, min(n + 1, len(self.coeffs)))
                    ),
                    ZERO,
                )
                values.append(-first * rest)
        return Jet(values)

    def divide_by_parameter(self):
        # This is an identity in Q(c), before either primitive specialization.
        assert not self or not self.coeffs[0]
        return Jet(self.coeffs[1:])


def certify(job, precision=11, limit=5):
    k, sign_w, color = job
    Jet.precision = precision
    relations = local_relations(
        Jet,
        k,
        1,
        sign_w,
        color,
        level=Jet.constant(Coefficient(C)),
        parameter=Jet([ZERO, ONE]),
    )
    check_cut(relations, Jet, job, limit)


def check_cut(relations, series, job, limit):
    cuts = local_central_sections(relations, series)
    center_precision = series.precision
    for start in range(4):
        series.precision = center_precision
        losses = []

        def project(n, h, rows, size):
            free, images, loss = saturated_projection(
                rows, size, DIMENSIONS[n][h], series
            )
            losses.append(loss)
            return free, images

        for n in quotient_layers(
            relations, cuts, start, limit, series.constant(1), project
        ):
            series.precision -= max(losses)
            losses.clear()
            assert series.precision > 0, (job, start, n)
        print(
            "Relative quartic center and cut",
            job,
            start,
            limit,
            series.precision,
            flush=True,
        )


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=2) as pool:
        list(pool.map(certify, product((1, -1), (1, -1), range(2))))
