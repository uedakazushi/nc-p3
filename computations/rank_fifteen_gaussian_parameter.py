"""Truncated Gaussian series with one rational-function parameter and recorded pivot norms."""

from flint import fmpz_mod_poly_ctx, fmpz_poly

from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_pole_slope_annulus import IntegralPoleJet


class SharedPrecision(type):
    @property
    def precision(cls):
        return cls.base.precision

    @precision.setter
    def precision(cls, digits):
        cls.base.precision = digits


def gaussian_parameter_jets(modulus):
    """Select units at (3, modulus) in Q(i)(m); divide series exactly over that field."""
    residual = fmpz_mod_poly_ctx(3)(list(modulus))
    assert residual.is_irreducible() and residual.degree() % 2 == 1

    class Base(IntegralPoleJet):
        pass

    Base.modulus = fmpz_poly(list(modulus))
    zero = Base([])

    def recorded_inverse(value):
        for polynomial in (value.coeffs[0], value.den):
            Recorded.record(tuple(int(c) % 3 for c in polynomial))
        return value.inverse()

    class Gaussian(metaclass=SharedPrecision):
        __slots__ = ("real", "imag")
        base = Base

        def __init__(self, real, imag=zero):
            self.real, self.imag = real, imag

        @classmethod
        def constant(cls, value):
            return cls(Base.constant(value))

        def __bool__(self):
            return bool(self.real) or bool(self.imag)

        def __neg__(self):
            return Gaussian(-self.real, -self.imag)

        def __add__(self, other):
            if not isinstance(other, Gaussian):
                other = Gaussian.constant(other)
            return Gaussian(self.real+other.real, self.imag+other.imag)

        __radd__ = __add__

        def __sub__(self, other):
            return self + -other

        def __rsub__(self, other):
            return -self + other

        def __mul__(self, other):
            if not isinstance(other, Gaussian):
                other = Gaussian.constant(other)
            if not self.imag:
                return Gaussian(self.real*other.real, self.real*other.imag)
            if not other.imag:
                return Gaussian(self.real*other.real, self.imag*other.real)
            return Gaussian(self.real*other.real-self.imag*other.imag,
                            self.real*other.imag+self.imag*other.real)

        __rmul__ = __mul__

        def unit(self):
            # An odd-degree extension of F3 still has -1 as a nonsquare.
            return self.real.unit() or self.imag.unit()

        def inverse(self):
            assert self.unit()
            if not self.imag:
                return Gaussian(recorded_inverse(self.real))
            if not self.real:
                return Gaussian(zero, -recorded_inverse(self.imag))
            reciprocal_norm = recorded_inverse(self.real*self.real+self.imag*self.imag)
            return Gaussian(self.real*reciprocal_norm, -self.imag*reciprocal_norm)

        def __truediv__(self, other):
            return self * other.inverse()

        def divide_by_parameter(self):
            return Gaussian(self.real.divide_by_parameter(), self.imag.divide_by_parameter())

    return Gaussian
