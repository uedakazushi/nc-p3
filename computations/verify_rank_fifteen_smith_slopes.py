"""Transfer the generic H-Smith exponents to alternative geometric residue points."""

from itertools import product

from flint import fmpz_mod_poly_ctx

from verify_rank_fifteen_integral_annulus import aligned_centers
from verify_rank_fifteen_integral_slope import trace_map
from verify_rank_fifteen_overlap_annulus import ONE, ZERO, TraceJet, poly
from verify_rank_fifteen_overlap_bound import overlap, prepend_tables, seed_layers
from verify_rank_fifteen_variable_zero_disk import trace_relations
from verify_rank_fourteen_all_ordinary import DIMENSIONS


def at_residue(coefficients):
    modulus = poly(list(coefficients))
    assert fmpz_mod_poly_ctx(3)(list(coefficients)).is_irreducible()
    orders = []

    class ResidueJet(TraceJet):
        """H-jets over F_3[z]/(modulus), with reduced polynomial coordinates."""

        record_pivot_orders = staticmethod(orders.append)

        def __init__(self, coefficients, den=None, *, reduce=True):
            coefficients = coefficients[:self.precision]
            if den is None or den == ONE:
                coefficients = [c % modulus for c in coefficients]
            else:
                gcd, inverse, _ = den.xgcd(modulus)
                assert gcd == ONE
                coefficients = [(c*inverse) % modulus for c in coefficients]
            while coefficients and not coefficients[-1]:
                coefficients.pop()
            self.coeffs, self.den = coefficients, ONE

        def unit(self):
            return bool(self.coeffs and self.coeffs[0])

        def inverse(self):
            gcd, first, _ = self.coeffs[0].xgcd(modulus)
            assert gcd == ONE
            values = [first]
            if len(self.coeffs) > 1:
                for k in range(1, self.precision):
                    total = sum((self.coeffs[j]*values[k-j]
                                 for j in range(1, min(k+1, len(self.coeffs)))), ZERO)
                    values.append((-total*first) % modulus)
            return ResidueJet(values)

    jet = ResidueJet

    def windows():
        return trace_relations(jet.constant(1), jet.constant(2), jet([poly([0, 1])]),
                               jet([ZERO, ONE]), jet.precision)

    relations, _, _, centers = aligned_centers(15, jet, windows)
    dims, left, right, words = {}, {}, {}, []
    unit = jet.constant(1)
    for p in range(4):
        jet.precision = 11
        words.append(seed_layers(p, 7, relations, centers, dims, right, jet,
                                 lambda value: jet([value.coeffs[0]])))
    for p in range(4):
        left.update(prepend_tables(p, 7, words, right, unit=unit))
    for n in range(8, 14):
        for p, h in product(range(4), repeat=2):
            dims[p, n, h] = overlap(n, p, h, dims, left, right, unit=unit)
        expected = (6,)*4 if n == 9 else DIMENSIONS[n]
        assert all(tuple(dims[p, n, h] for h in range(4)) == expected for p in range(4))
    assert len(orders) == 200
    return tuple(orders)


# Minimal polynomials for D and for a selected z above it, in ascending degree.
CHARTS = (
    ((2, 1), (2, 0, 1, 2, 2, 0, 1)),
    ((1, 1), (1, 1, 1, 1, 1)),
    ((2, 2, 1), (2, 0, 0, 2, 1)),
)


def fiber_polynomial(z, coefficients):
    p, q = trace_map(z)
    degree = len(coefficients)-1
    return sum((c*p**j*q**(degree-j) for j, c in enumerate(coefficients)), 0*z)


def verify():
    ring = fmpz_mod_poly_ctx(3)
    reference = at_residue((2, 1, 1))
    for g, m in CHARTS:
        modulus = ring(list(m))
        assert ring(list(g)).is_irreducible()
        assert not fiber_polynomial(ring.gen(), g) % modulus
        assert at_residue(m) == reference
        print("All 200 Smith profiles and the residue overlap verified for", g, flush=True)


if __name__ == "__main__":
    verify()
