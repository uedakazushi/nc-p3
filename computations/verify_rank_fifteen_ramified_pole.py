"""Resolve the paired square roots where x has residue one on the pole chart."""

from flint import ctx, fmpq_poly, fmpz_mod_poly_ctx, fmpz_poly
import sympy as sp

from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_pole_fibers import cuts, tables, transform
from verify_rank_fifteen_pole_slope_annulus import (
    IntegralPoleJet, ONE, TraceJet, initial_tables, pattern, poly, tensor_relations,
)

from verify_rank_fourteen_all_ordinary import gaussian_jets

LOSSES = (1, 0, 1, 1, 2, 0, 1, 0)
MODULUS = poly([2, 1, 1])


class RamifiedPole(TraceJet):
    def unit(self):
        assert self.den % MODULUS
        if not self.coeffs or not self.coeffs[0] % MODULUS:
            return False
        for coefficient in (self.coeffs[0], self.den):
            Recorded.record(tuple(map(int, coefficient)))
        return True


class ExactRamifiedPole(IntegralPoleJet):
    modulus = fmpz_poly([2, 1, 1])


# Polynomials for the tangent value J and a chosen z above each residual orbit.
EXCEPTIONS = (
    ((1, 1), (2, 0, 1, 0, 1)),
    ((1, 0, 1, 1, 1), (2, 1, 0, 0, 1)),
    ((1, 0, 1), (1, 1, 1, 1, 1)),
    ((2, 1, 1), (1, 0, 1, 2, 1)),
    ((2, 0, 1, 1), (1, 2, 0, 2, 0, 2, 1)),
)


class Supplemental(TraceJet):
    moduli = tuple(poly(list(m)) for _, m in EXCEPTIONS)

    def unit(self):
        assert all(self.den % m for m in self.moduli)
        return bool(self.coeffs) and all(self.coeffs[0] % m for m in self.moduli)


def relations(jet, half, h, parameter):
    one = jet.constant(1)
    denominator = (one-2*h*h).inverse()
    inverse = parameter.inverse()
    along = h*(parameter+inverse)*half
    across = h*(parameter-inverse)*half

    def root(square, initial):
        value = initial
        for _ in range(jet.precision.bit_length()):
            value = (value+square*value.inverse())*half
        assert not (value*value-square)
        return value

    e = root(one-along*along, one)
    f = root(one-along*along-h*h, one)
    n, b = (2*h*denominator*(e+sign*along) for sign in (-1, 1))
    r2 = 8*h*h*denominator*denominator
    t = (one-2*along*along)*half
    assert not (n*b-r2*t) and not (n*n+b*b-r2)
    p, u = (root(4*denominator+sign*2*t*r2, jet.constant(2)) for sign in (1, -1))
    q, v = 4*h*denominator*f, 4*h*denominator*across
    gamma = pattern(one, 2*h*h*denominator, n*half, b*half, (u-p)*half*half,
                    (u+p)*half*half, (v-q)*half*half, (v+q)*half*half)
    return tensor_relations(gamma, one*0)


def compute(jet, precision, *, exact=False, last=13):
    orders = []

    class Profile(jet):
        record_pivot_orders = staticmethod(orders.append)

    products, coordinates = tables(Profile, precision, exact=exact, build=relations,
                                  uncut_losses=(0, 2, 2, 2, 2), fifth_loss=1)
    cuts(Profile, products, last=last, losses=LOSSES)
    assert Profile.precision > 0 and len(orders) == 180 + 16*(last-5)
    return coordinates, tuple(orders)


def verify_tangents(divisor):
    ring = fmpz_mod_poly_ctx(3)
    z = ring.gen()
    num, den = (z**4+1)**2, 4*z**4
    core = ring(list(map(int, divisor)))
    # Zero and infinity both lie over the pole of the rational map.
    while not core[0]:
        core = core // z
    while core.degree() > 0:
        degree = core.degree()
        core = core.gcd(core(-z)).monic()
        assert all(not core[j] for j in range(1, core.degree()+1, 2))
        rotated = ring([int(core[j])*(-1)**(j//2) for j in range(core.degree()+1)])
        core = core.gcd(rotated).monic()
        core = core.gcd(core.reverse(degree=core.degree())).monic()
        if core.degree() == degree:
            break
    # J=1 is a branch fiber with double roots.
    assert num-den == (z**4-1)**2
    expected = z**4-1
    for g, m in EXCEPTIONS:
        polynomial, modulus = ring(list(g)), ring(list(m))
        fiber = transform(polynomial, num, den)
        assert polynomial.is_irreducible() and modulus.is_irreducible()
        assert not fiber % modulus and den % modulus
        expected *= fiber
    assert core.degree() == 100 and core == expected.monic()


def collision():
    a, ell = sp.symbols("a ell")
    t = 1-2*a*ell
    scaled_r = 2*a/(1-2*a)**2
    x = (1+scaled_r)*t**2
    root_y = scaled_r*t**3/2
    factors = (
        -4*a*(ell-1)*(a*ell-1)*(4*a**2*ell-1),
        4*a*ell*(a*ell+a-1)*(4*a**2*ell-4*a+1),
    )
    for sign, factor in zip((1, -1), factors):
        assert sp.cancel((x-1+sign*2*root_y)*(1-2*a)**2-factor) == 0
    ctx.cap = 25
    jet = gaussian_jets(-1, 20, residue_prime=3)
    h = jet(fmpq_poly([0, 1]))
    i = jet(fmpq_poly([]), fmpq_poly([1]))
    products, _ = initial_tables(
        jet, 20, jet(fmpq_poly([1])/2), h, i, build=relations,
        uncut_losses=(0, 2, 2, 3, 2), center_losses=(0, 1, 0, 0), fourth_loss=1, fifth_loss=1,
    )
    cuts(jet, products, losses=(1, 0, 1, 1, 2, 0, 3, 0))
    assert jet.precision == 1


def verify():
    Recorded.divisor = ONE
    Recorded.record.cache_clear()
    coordinates, reference = compute(RamifiedPole, 17)
    print("Tangent pivot divisor degree", Recorded.divisor.degree(), flush=True)
    assert Recorded.divisor.degree() == 897
    verify_tangents(Recorded.divisor)
    assert compute(Supplemental, 17) == (coordinates, reference)
    exact_coordinates, orders = compute(ExactRamifiedPole, 12, exact=True, last=7)
    assert exact_coordinates == coordinates and orders == reference[:len(orders)]
    collision()
    print("Verified the six-digit tangent restart and the unrestricted-depth collision tube.")


if __name__ == "__main__":
    verify()
