"""An unrestricted-depth certificate on x=1 over Q(zeta_8), localized above 3."""

from functools import cache

import sympy as sp
from flint import fmpq, fmpq_poly

from verify_rank_fifteen_generic_levels import CommonDenominatorJet
from verify_rank_fifteen_pole_fibers import cuts
from verify_rank_fifteen_pole_slope_annulus import initial_tables
from verify_rank_fifteen_ramified_pole import LOSSES, relations

ONE = fmpq_poly([1])
ZERO = fmpq_poly([])
Z = fmpq_poly([0, 1])
MODULUS = Z**4+1


@cache
def hensel_coefficient(modulus):
    """Lift c=1 mod 3 with c^2=-2 to the requested power of three."""
    c = 1
    while (c*c+2) % modulus:
        c = (c-(c*c+2)*pow(2*c, -1, modulus)) % modulus
    return c


def local_unit(value):
    denominator = int(value.denom())
    power = 1
    while denominator % 3 == 0:
        denominator //= 3
        power *= 3
    modulus = 3*power
    c = hensel_coefficient(modulus)
    # Modulo Z^2+cZ-1, Z^2=1-cZ and Z^3=-c-Z.
    g = value.numer()
    u = int(g[0])+int(g[2])-c*int(g[3])
    v = int(g[1])-c*int(g[2])-int(g[3])
    assert u % power == v % power == 0
    return bool(u % modulus or v % modulus)


class CyclotomicJet(CommonDenominatorJet):
    poly_one, poly_zero = ONE, ZERO

    def __init__(self, coeffs, den=None, *, reduce=True):
        assert den is None or den == ONE
        values = [c % MODULUS for c in coeffs[:self.precision]]
        while values and not values[-1]:
            values.pop()
        self.coeffs, self.den = values, ONE

    @classmethod
    def constant(cls, value):
        return cls([fmpq_poly([value])])

    def unit(self):
        return bool(self.coeffs) and local_unit(self.coeffs[0])

    def inverse(self):
        assert self.unit()
        gcd, first, _ = self.coeffs[0].xgcd(MODULUS)
        assert gcd == ONE
        values = [first]
        if len(self.coeffs) > 1:
            for n in range(1, self.precision):
                total = sum((self.coeffs[j]*values[n-j]
                             for j in range(1, min(n+1, len(self.coeffs)))), ZERO)
                values.append((-first*total) % MODULUS)
        return type(self)(values)


def verify():
    t = sp.symbols("t", nonzero=True)
    p, q = (t+1/t)/2, (t-1/t)/(3*sp.sqrt(3))
    assert sp.cancel(p**2-27*q**2/4-1) == 0
    assert sp.cancel(p+3*sp.sqrt(3)*q/2-t) == 0
    jet = CyclotomicJet
    jet.precision = 17
    h, zeta = jet([ZERO, ONE]), jet([Z])
    one, half = jet.constant(1), jet.constant(fmpq(1, 2))
    d = (one-2*h*h).inverse()

    def root(square, initial):
        value = initial
        for _ in range(jet.precision.bit_length()):
            value = (value+square*value.inverse())*half
        assert not (value*value-square)
        return value

    g = root(one+2*h*h*d*d, one)
    l2 = d*d*(g*(g+one)).inverse()
    ell = root(l2, (zeta+zeta.inverse())*half)
    across = root(l2-one, (zeta-zeta.inverse())*half)
    assert not (g*g*(one-2*h*h*ell*ell)*(one-2*h*h*ell*ell)-one)
    products, _ = initial_tables(jet, 17, half, h, ell+across, build=relations,
                                 uncut_losses=(0, 2, 2, 2, 2), fifth_loss=1)
    cuts(jet, products, losses=LOSSES)
    assert jet.precision == 1
    print("Verified the vertical conic and its seventeen-digit cyclotomic certificate.")


if __name__ == "__main__":
    verify()
