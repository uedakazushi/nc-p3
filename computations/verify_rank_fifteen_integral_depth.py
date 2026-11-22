"""Descend the three-residue-class certificate at unrestricted depth in characteristic three."""

from itertools import product

from flint import fmpz_mod_poly_ctx

from verify_rank_fifteen_integral_annulus import Recorded, windows as trace_windows
from verify_rank_fifteen_integral_slope import invariant_core
from verify_rank_fifteen_overlap_annulus import ONE, TraceJet, poly
from verify_rank_fifteen_overlap_bound import prepend_tables, seed_layers
from verify_rank_fifteen_smith_slopes import fiber_polynomial
from verify_rank_fifteen_two_sided_cut import overlap_rows, store_projection
from verify_rank_fifteen_variable_zero_disk import center_data
from verify_rank_fourteen_all_ordinary import DIMENSIONS, saturated_projection


def residues(polynomial):
    return tuple(int(c) % 3 for c in polynomial)


class Symbolic(TraceJet):
    def unit(self):
        if not super().unit():
            return False
        Recorded.record(residues(self.coeffs[0]))
        Recorded.record(residues(self.den))
        return True


def constant_term(value):
    return Recorded([poly(list(residues(value.coeffs[0])))], poly(list(residues(value.den))))


def compute_chart():
    Recorded.divisor = ONE
    Recorded.record.cache_clear()
    windows, _, _, centers = center_data(15, Symbolic, lambda: trace_windows(Symbolic))
    Recorded.precision = 1
    dims, right, left, words = {}, {}, {}, []
    for p in range(4):
        Symbolic.precision = 11
        words.append(seed_layers(p, 7, windows, centers, dims, right, Symbolic, constant_term,
                                 keep=lambda value: value and value.coeffs[0]))
        assert Symbolic.precision == 1
        print("Recorded length-seven pivots at start", p, flush=True)
    for p in range(4):
        left.update(prepend_tables(p, 7, words, right, unit=Recorded.constant(1)))
    for n in range(8, 14):
        for p, h in product(range(4), repeat=2):
            ends, beginnings, width, rows = overlap_rows(n, p, h, dims, left, right)
            target = 6 if n == 9 else DIMENSIONS[n][h]
            free, projection, loss = saturated_projection(rows, width, target, Recorded)
            assert loss == 0
            store_projection(n, p, h, dims, left, right, ends, beginnings, free, projection)
        print("Recorded residue overlaps at length", n, flush=True)
    return Recorded.divisor


def slope_core(divisor):
    ring = fmpz_mod_poly_ctx(3)
    z = ring.gen()
    f = ring(list(map(int, divisor)))
    return invariant_core(f, z, ring.one(), lambda a, b: a.gcd(b).monic(), lambda a: a.monic())


def verify():
    divisor = compute_chart()
    assert divisor.degree() == 12017
    core = slope_core(divisor)
    ring = fmpz_mod_poly_ctx(3)
    expected = ring.one()
    for g in ((1, 1), (2, 1), (2, 2, 1)):
        expected *= fiber_polynomial(ring.gen(), g)
    assert core.degree() == 96 and core == expected.monic()
    print("All primary failures lie over the two prime-field units and one conjugate pair.")


if __name__ == "__main__":
    verify()
