"""Four-central-quartic cuts on the integral-slope three-adic annulus."""

from functools import lru_cache

from flint import fmpz_mod_poly_ctx

from verify_rank_fifteen_integral_slope import RecordedJet as IntegralJet, invariant_core, trace_map
from verify_rank_fifteen_overlap_annulus import ONE, ZERO, TraceJet, poly
from verify_rank_fifteen_sixth_center_minors import verify as verify_sixth
from verify_rank_fifteen_two_sided_cut import verify as verify_overlap
from verify_klein_center_lift import central_equations
from verify_rank_fifteen_variable_zero_disk import FIXED_CENTER_WORDS, center_data, trace_relations
from verify_rank_fourteen_all_ordinary import saturated_projection

# (Minimal polynomial of the residue slope, one minimal polynomial for its trace parameter.)
CHARTS = (((1, 1), (1, 1, 1, 1, 1)),
          ((2, 2, 1), (2, 0, 0, 2, 1)),
          ((2, 2, 2, 1), (2, 2, 0, 1)))


class Recorded(TraceJet):
    divisor = ONE

    @classmethod
    @lru_cache(maxsize=12000)
    def record(cls, coefficients):
        f = poly(list(coefficients))
        while f.degree() > 0:
            common = f.gcd(cls.divisor)
            if common.degree() == 0:
                for factor, _ in f.factor_squarefree()[1]:
                    cls.divisor *= factor
                return
            f = f // common

    def unit(self):
        if not super().unit():
            return False
        self.record(tuple(map(int, self.coeffs[0])))
        self.record(tuple(map(int, self.den)))
        return True


class Supplemental(TraceJet):
    moduli = tuple(poly(list(m)) for _, m in CHARTS)

    def unit(self):
        assert all(self.den % m for m in self.moduli)
        return bool(self) and all(self.coeffs[0] % m for m in self.moduli)


def windows(jet):
    return trace_relations(jet.constant(1), jet.constant(2), jet([poly([0, 1])]),
                           jet([ZERO, ONE]), jet.precision)


def aligned_centers(precision, jet, build_relations):
    windows, basis, image, _ = center_data(precision, jet, build_relations, fixed_sections=False)
    center_basis = basis.copy()
    centers = [[] for _ in range(4)]
    for h, dimension in enumerate((3, 2)):
        words = [tuple(map(int, label)) for label in FIXED_CENTER_WORDS[h]]
        _, _, loss = saturated_projection([image(3, w).copy() for w in words], len(words), 0, jet)
        assert loss == 0
        center_basis[3, 4, h] = words
        rows, width = central_equations(center_basis, image, h)

        def choose(rows):
            return next(((i, k) for i in sorted(range(len(rows)), key=lambda j: len(rows[j]))
                         for k in sorted(rows[i]) if k < width-dimension and rows[i][k].unit()), None)

        free, projection, loss = saturated_projection(rows, width, dimension, jet, select_pivot=choose)
        assert loss == 0 and free == list(range(width-dimension, width))
        for j in range(1 if h == 0 else 0, dimension):
            offset = 0
            for p in range(4):
                paths = center_basis[p, 4, h]
                vector = {w: projection[offset+i][j] for i, w in enumerate(paths)
                          if j in projection[offset+i]}
                centers[p].append((h, vector))
                offset += len(paths)
    return windows, basis, image, centers


def verify():
    verify_sixth(IntegralJet)
    assert IntegralJet.divisor.degree() == 205
    verify_overlap(Recorded, lambda: windows(Recorded))
    assert Recorded.divisor.degree() == 9835
    ring = fmpz_mod_poly_ctx(3)
    z, one = ring.gen(), ring.one()
    sixth = ring(list(IntegralJet.divisor))
    high = ring(list(map(int, Recorded.divisor)))
    core = invariant_core(sixth*high, z, one, lambda a, b: a.gcd(b).monic(), lambda f: f.monic())
    numerator, denominator = trace_map(z)
    exceptional = one
    for coefficients, m in CHARTS:
        g, modulus = ring(list(coefficients)), ring(list(m))
        fiber = sum((g[j]*numerator**j*denominator**(g.degree()-j)
                     for j in range(g.degree()+1)), one*0)
        assert g.is_irreducible() and modulus.is_irreducible()
        assert sixth % modulus and denominator % modulus and not fiber % modulus
        exceptional *= fiber
    assert core.degree() == 144 and core == exceptional.monic()
    verify_overlap(Supplemental, lambda: windows(Supplemental), center_builder=aligned_centers)
    print("All integral slopes, including every geometric root of all three residual polynomials.")


if __name__ == "__main__":
    verify()
