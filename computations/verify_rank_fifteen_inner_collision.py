"""Exact mixed-characteristic charts at z=1 inside the ramified pole collision."""

from flint import fmpz_poly, nmod_poly

from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_pole_fibers import RecordingExact, cuts, tables
from verify_rank_fifteen_ramified_pole import relations

T = nmod_poly([0, 1], 3)
EXTRA = ((2, 1, 1), (1, 0, 1), (1, 0, 2, 1))


class Supplemental(RecordingExact):
    moduli = tuple(fmpz_poly(list(coefficients)) for coefficients in EXTRA)

    def unit(self):
        def good(value):
            return all(any(int(c) % 3 for c in value % modulus) for modulus in self.moduli)
        assert good(self.den)
        if not self.coeffs or not good(self.coeffs[0]):
            return False
        for coefficient in (self.coeffs[0], self.den):
            Recorded.record(tuple(int(c) % 3 for c in coefficient))
        return True


def divisor(jet, power):
    Recorded.divisor = 1+0*T
    Recorded.record.cache_clear()
    precision = 16 if power == 2 else 15

    def collision(field, half, h, m):
        displacement = m*h*h if power == 2 else m*h
        return relations(field, half, h, field.constant(1)+displacement)

    products, _ = tables(
        jet, precision, exact=True, build=collision, uncut_losses=(0, 2, 2, 3, 2),
        center_losses=(int(power == 2), 0, 0, 0), fifth_loss=0,
    )
    cuts(jet, products, losses=(2, 0, 1, 0, 1, 0, 1, 0))
    assert jet.precision == 1
    return Recorded.divisor


def verify():
    class Deep(RecordingExact):
        modulus = fmpz_poly([0, 1])

    class Middle(RecordingExact):
        modulus = fmpz_poly([-1, 1])

    assert divisor(Deep, 2) == T-1
    baseline = divisor(Middle, 1)
    assert baseline.degree() == 20
    conjugate = nmod_poly([(-1 if i % 2 else 1)*int(c) for i, c in enumerate(baseline)], 3)
    core = baseline.gcd(conjugate)
    supplement = divisor(Supplemental, 1)
    covered = T
    for coefficients in EXTRA:
        factor = nmod_poly(list(coefficients), 3)
        reflected = nmod_poly([(-1 if i % 2 else 1)*int(c) for i, c in enumerate(factor)], 3)
        assert supplement.gcd(factor).degree() == 0
        covered *= factor*reflected
    assert covered % core == 0
    print("Verified z=1 collision for nu(J-1)>=nu(y) and 2nu(J-1)=nu(y).")


if __name__ == "__main__":
    verify()
