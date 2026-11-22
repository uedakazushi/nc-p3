"""Descend two-sided pole-chart pivots along the degree-eight conic map."""

from itertools import product

from flint import fmpz_mod_poly_ctx, fmpz_poly

from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_pole_slope_annulus import (
    IntegralPoleJet, ONE, PoleJet, TraceJet, ZERO, initial_tables, overlap, poly,
)

LOSSES = (2, 1, 1, 1, 1, 0, 1, 0)
# Each entry pairs a residual x-polynomial with a chosen lift in the s-coordinate.
EXCEPTIONS = (((2, 2, 2, 1), (2, 1, 1, 1)),
              ((2, 0, 1, 1), (2, 0, 1, 2, 2, 0, 1)))


class RecordingPole(PoleJet):
    def unit(self):
        if not super().unit():
            return False
        for coefficient in (self.coeffs[0], self.den):
            Recorded.record(tuple(map(int, coefficient)))
        return True


class RecordingExact(IntegralPoleJet):
    def unit(self):
        if not super().unit():
            return False
        for coefficient in (self.coeffs[0], self.den):
            Recorded.record(tuple(int(c) % 3 for c in coefficient))
        return True


class SupplementalPole(TraceJet):
    moduli = tuple(poly(list(m)) for _, m in EXCEPTIONS)

    def unit(self):
        assert all(self.den % m for m in self.moduli)
        return bool(self.coeffs) and all(self.coeffs[0] % m for m in self.moduli)


def tables(jet, precision, *, exact=False, **options):
    jet.precision = precision
    if exact:
        r = jet([fmpz_poly([]), fmpz_poly([1])])
        s = jet([fmpz_poly([0, 1])])
        half = jet([fmpz_poly([1])], fmpz_poly([2]))
    else:
        r, s, half = jet([ZERO, ONE]), jet([poly([0, 1])]), jet.constant(2)
    return initial_tables(jet, precision, half, r, s, **options)


def cuts(jet, products, last=13, losses=LOSSES):
    dims, left, right = products
    for n in range(6, last+1):
        used = [overlap(n, p, h, dims, left, right, jet) for p, h in product(range(4), repeat=2)]
        assert max(used) == losses[n-6]
        jet.precision -= max(used)
        assert jet.precision > 0


def slope(s):
    return 16*s**2*(1-s**2)**2, (1+s**2)**4


def transform(f, a, b):
    value = 0*a + f.leading_coefficient()
    denominator = 0*a + 1
    for i in reversed(range(f.degree())):
        denominator *= b
        value = value*a + f[i]*denominator
    return value


def check_fibers(divisor):
    ring = fmpz_mod_poly_ctx(3)
    s = ring.gen()
    num, den = slope(s)
    core = ring(list(map(int, divisor)))
    for branch in (num, den, num-den):
        while (factor := core.gcd(branch)).degree() > 0:
            core = core // factor
    while core.degree() > 0:
        prior = core.degree()
        for a, b in ((-s, ring.one()), (ring.one(), s), (1+s, 1-s)):
            assert transform(num, a, b)*b**2*den == transform(den, a, b)*num
            core = core.gcd(transform(core, a, b)).monic()
        if core.degree() == prior:
            break
    expected = ring.one()
    for g, m in EXCEPTIONS:
        xpoly, modulus = ring(list(g)), ring(list(m))
        fiber = transform(xpoly, num, den)
        assert xpoly.is_irreducible() and modulus.is_irreducible()
        assert not fiber % modulus and den % modulus
        expected *= fiber
    assert core.degree() == 48 and core == expected.monic()


def verify():
    Recorded.divisor = ONE
    Recorded.record.cache_clear()
    exact_tables, exact_coordinates = tables(RecordingExact, 12, exact=True)
    cuts(RecordingExact, exact_tables, last=6)
    assert RecordingExact.precision == 1
    s = poly([0, 1])
    num, den = slope(s)
    assert Recorded.divisor == s**9-s and not num*den*(num-den) % Recorded.divisor
    Recorded.divisor = ONE
    Recorded.record.cache_clear()
    products, coordinates = tables(RecordingPole, 17)
    assert coordinates == exact_coordinates
    cuts(RecordingPole, products)
    assert Recorded.divisor.degree() == 121
    check_fibers(Recorded.divisor)
    products, coordinates = tables(SupplementalPole, 17)
    assert coordinates == exact_coordinates
    cuts(SupplementalPole, products)
    print("Verified the two cubic residue fibers and the uniform six-digit restart.")


if __name__ == "__main__":
    verify()
