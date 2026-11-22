"""Descend the sixth central certificate through the degree-24 trace map."""

from functools import lru_cache

from flint import fmpz_mod_poly_ctx, fmpz_poly

from verify_rank_fifteen_sixth_center_minors import TraceJet, verify as verify_minors

ONE = fmpz_poly([1])
Z = fmpz_poly([0, 1])
PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 53, 59,
          71, 79, 89, 131, 173, 179, 191, 211, 233, 271, 313, 317, 349, 373, 379,
          433, 503, 911, 1373)
EXCEPTIONS = {5: 0, 7: 0, 13: 0, 19: 17, 29: 3, 37: 0}


# One irreducible residue factor above each fiber retained by the first chart.
RESIDUE_CHARTS = ((5, (2, 1, 1)), (7, (6, 1, 1)), (13, (11, 1)),
                  (29, (27, 1)), (37, (35, 1)))


class FiberJet(TraceJet):
    rings = tuple((fmpz_mod_poly_ctx(p), fmpz_mod_poly_ctx(p)(list(m)))
                  for p, m in RESIDUE_CHARTS)

    @classmethod
    @lru_cache(maxsize=8192)
    def units(cls, coefficients):
        return all(bool(ring(list(coefficients)) % modulus) for ring, modulus in cls.rings)

    def unit(self):
        assert self.units(tuple(self.den))
        return bool(self.coeffs) and self.units(tuple(self.coeffs[0]))


class Fiber19Jet(FiberJet):
    rings = ((fmpz_mod_poly_ctx(19), fmpz_mod_poly_ctx(19)([17, 1])),)


def verify_remaining_fibers():
    for jet, primes in ((FiberJet, tuple(p for p, _ in RESIDUE_CHARTS)), (Fiber19Jet, (19,))):
        verify_fibers(jet, primes)


def verify_fibers(jet, primes):
    for p, (ring, modulus) in zip(primes, jet.rings):
        z = ring.gen()
        numerator, denominator = trace_map(z)
        assert modulus.is_irreducible() and denominator % modulus
        assert not (numerator-EXCEPTIONS[p]*denominator) % modulus
    verify_minors(jet, fixed_sections=False)


def trace_map(z):
    return (z**8+14*z**4+1)**3, 4*z**4*(z**4-1)**4


def verify_symmetries():
    # Multiplication by i fixes z^4. Clear the denominators for the other generator.
    p, q = trace_map(Z)
    a, b = Z+1, Z-1
    transformed_p = (a**8+14*a**4*b**4+b**8)**3
    transformed_q = 4*a**4*b**4*(a**4-b**4)**4
    assert p*transformed_q == q*transformed_p


class RecordedJet(TraceJet):
    divisor = ONE

    @classmethod
    def record(cls, poly):
        content = abs(int(poly.content()))
        assert content > 0 and content & (content-1) == 0
        f = poly // content
        if f.leading_coefficient() < 0:
            f = -f
        while f.degree() > 0:
            common = cls.divisor.gcd(f)
            if common == ONE:
                cls.divisor *= f // f.gcd(f.derivative())
                break
            f = f // common

    def unit(self):
        ok = super().unit()
        if ok:
            self.record(self.coeffs[0])
            self.record(self.den)
        return ok


def invariant_core(f, z, one, gcd, normalize):
    pole = z*(z**4-1)
    while (d := f.gcd(pole)).degree() > 0:
        f = f // d
    f = normalize(f)
    while f.degree() > 0:
        degree = f.degree()
        f = gcd(f, f(-z))
        assert all(not f[j] for j in range(1, f.degree()+1, 2))
        twist = sum((f[j]*(-1)**(j//2)*z**j for j in range(0, f.degree()+1, 2)), z*0)
        f = gcd(f, twist)

        # Homogenized substitution z -> (z+1)/(z-1).
        transformed = one*f.leading_coefficient()
        denominator = one
        for j in reversed(range(f.degree())):
            denominator *= z-1
            transformed = transformed*(z+1) + f[j]*denominator
        f = gcd(f, transformed)
        if f.degree() == degree:
            break
    return f


def supported(integer):
    n = abs(int(integer))
    assert n
    for p in PRIMES:
        while n % p == 0:
            n //= p
    assert n == 1


def integer_gcd(a, b):
    d = a.gcd(b)
    supported(a.leading_coefficient())
    supported(b.leading_coefficient())
    supported(d.content())
    supported((a//d).resultant(b//d))
    return d // d.content()


def verify():
    verify_symmetries()
    RecordedJet.divisor = ONE
    verify_minors(RecordedJet)
    f = RecordedJet.divisor
    print("Pivot divisor degree", f.degree(), flush=True)
    assert f.degree() == 205
    assert invariant_core(f, Z, ONE, integer_gcd, lambda g: g) == ONE
    for p in PRIMES[1:]:
        ring = fmpz_mod_poly_ctx(p)
        z, one = ring.gen(), ring.one()
        core = invariant_core(ring(list(f)), z, one, lambda a, b: a.gcd(b).monic(), lambda g: g.monic())
        if p in EXCEPTIONS:
            assert core.degree() == (24 if p in (13, 19, 29) else 8)
            numerator, denominator = trace_map(z)
            fiber = numerator-EXCEPTIONS[p]*denominator
            while (d := core.gcd(fiber)).degree() > 0:
                core = core // d
        assert core.degree() == 0
    verify_remaining_fibers()
    print("All integral slopes covered in residue characteristic zero or odd.")


if __name__ == "__main__":
    verify()
