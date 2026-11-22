"""Mixed-characteristic four-cut when 3 vanishes modulo H^15."""

from flint import nmod_poly

from verify_rank_fifteen_generic_levels import CommonDenominatorJet
from verify_rank_fifteen_two_sided_cut import verify
from verify_rank_fifteen_variable_zero_disk import trace_relations


def poly(cs):
    return nmod_poly(cs, 3)


ONE = poly([1])
ZERO = poly([])
MODULUS = poly([2, 1, 1])


class Fraction:
    """A rational function over F_3 with a polynomial numerator and denominator."""

    __slots__ = ('num', 'den')

    def __init__(self, num, den=ONE):
        if not den:
            raise ZeroDivisionError
        if not num:
            den = ONE
        elif den != ONE:
            common = num.gcd(den)
            num, den = num//common, den//common
        self.num, self.den = num, den

    @classmethod
    def constant(cls, value):
        return cls(poly([int(value) % 3]))

    def __bool__(self):
        return bool(self.num)

    def __neg__(self):
        return Fraction(-self.num, self.den)

    def __add__(self, other):
        if not self:
            return other
        if not other:
            return self
        if self.den == other.den:
            return Fraction(self.num+other.num, self.den)
        common = self.den.gcd(other.den)
        a, b = self.den//common, other.den//common
        return Fraction(self.num*b+other.num*a, self.den*b)

    def __mul__(self, other):
        return Fraction(self.num*other.num, self.den*other.den)

    def inverse(self):
        return Fraction(self.den, self.num)


class TraceJet(CommonDenominatorJet):
    poly_one = ONE
    poly_zero = ZERO
    coefficient = Fraction
    coefficient_zero = Fraction(ZERO)
    signed = False

    def unit(self):
        assert self.den % MODULUS
        return bool(self) and bool(self.coeffs[0] % MODULUS)


def relations():
    return trace_relations(TraceJet.constant(1), TraceJet.constant(2), TraceJet([poly([0, 1])]),
                           TraceJet([ZERO, ONE]), TraceJet.precision)


if __name__ == '__main__':
    verify(TraceJet, relations)
