"""Exact rational-function arithmetic using python-flint polynomial gcd.

Each element has a reduced numerator and monic denominator. No floating point,
finite-field evaluations, or interpolation are used. Independent field contexts
can coexist; the default field is QQ(u,v).
"""
from functools import lru_cache
from flint import fmpq_mpoly_ctx, fmpq, nmod_mpoly_ctx
import sympy as s


def rational_function_field(names=('u', 'v'), *, prime=None):
    if prime is None:
        ctx = fmpq_mpoly_ctx.get(tuple(names))
    else:
        if not s.isprime(prime):
            raise ValueError("The coefficient modulus must be prime")
        ctx = nmod_mpoly_ctx.get(tuple(names), prime)

    def scalar(value):
        value = fmpq(value)
        if prime is None:
            return value
        denominator = int(value.denominator) % prime
        if denominator == 0:
            raise ZeroDivisionError
        return int(value.numerator) * pow(denominator, -1, prime) % prime

    variables = s.symbols(tuple(names))
    unit = ctx.constant(1)

    class Element:
        __slots__ = ('n', 'd')

        def __init__(self, numerator, denominator=unit):
            if not denominator:
                raise ZeroDivisionError
            if not numerator:
                self.n, self.d = ctx.constant(0), unit
                return
            if denominator != unit:
                common = numerator.gcd(denominator)
                numerator, denominator = numerator // common, denominator // common
                lead = denominator.leading_coefficient()
                if lead != 1:
                    numerator, denominator = numerator / lead, denominator / lead
            self.n, self.d = numerator, denominator

        @staticmethod
        def coerce(other):
            if isinstance(other, Element):
                return other
            return Element(ctx.constant(scalar(other)))

        def __bool__(self):
            return bool(self.n)

        def __eq__(self, other):
            other = self.coerce(other)
            return self.n == other.n and self.d == other.d

        def __add__(self, other):
            other = self.coerce(other)
            if not other:
                return self
            if not self:
                return other
            return Element(self.n * other.d + other.n * self.d, self.d * other.d)

        __radd__ = __add__

        def __neg__(self):
            return Element(-self.n, self.d)

        def __sub__(self, other):
            return self + -self.coerce(other)

        def __rsub__(self, other):
            return self.coerce(other) + -self

        def __mul__(self, other):
            other = self.coerce(other)
            if not self or not other:
                return zero
            return Element(self.n * other.n, self.d * other.d)

        __rmul__ = __mul__

        def __truediv__(self, other):
            other = self.coerce(other)
            return Element(self.n * other.d, self.d * other.n)

        def __rtruediv__(self, other):
            return self.coerce(other) / self

        def __pow__(self, exponent):
            if exponent < 0:
                return Element(self.d**(-exponent), self.n**(-exponent))
            return Element(self.n**exponent, self.d**exponent)

        def __repr__(self):
            return f'({self.n})/({self.d})'

        def as_expr(self):
            def polynomial(value):
                return sum((s.Rational(str(coefficient)) * s.prod(
                    variable ** power for variable, power in zip(variables, monomial))
                    for monomial, coefficient in value.to_dict().items()), s.S.Zero)
            return polynomial(self.n) / polynomial(self.d)

    zero, one = Element(ctx.constant(0)), Element(unit)

    class Field:
        @staticmethod
        @lru_cache(None)
        def from_sympy(expression):
            numerator, denominator = s.fraction(s.cancel(expression))
            def polynomial(value):
                return ctx.from_dict({monomial: scalar(fmpq(int(coefficient.p), int(coefficient.q)))
                    for monomial, coefficient in s.Poly(value, *variables).terms()})
            return Element(polynomial(numerator), polynomial(denominator))

    result = Field()
    result.zero, result.one = zero, one
    return result


K = rational_function_field()
zero, one = K.zero, K.one
