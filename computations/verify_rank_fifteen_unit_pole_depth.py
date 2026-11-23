"""Exact mixed-characteristic lattices at arbitrary depth on both unit pole charts."""

import argparse

from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_pole_fibers import (
    ONE, RecordingExact, RecordingPole, SupplementalPole, check_fibers, cuts, tables,
)
from verify_rank_fifteen_ramified_pole import (
    ExactRamifiedPole, RamifiedPole, Supplemental, compute, verify_tangents,
)


def reset_divisor():
    Recorded.divisor = ONE
    Recorded.record.cache_clear()


def conic_profile(jet, *, exact=False):
    orders = []

    class Profile(jet):
        record_pivot_orders = staticmethod(orders.append)

    products, coordinates = tables(Profile, 17, exact=exact)
    cuts(Profile, products)
    assert Profile.precision == 1 and len(orders) == 308
    return coordinates, orders


def nonbranch():
    reset_divisor()
    reference = conic_profile(RecordingExact, exact=True)
    exact_divisor = Recorded.divisor
    print("Exact nonbranch pole lattices through thirteen.", flush=True)
    reset_divisor()
    assert conic_profile(RecordingPole) == reference
    assert Recorded.divisor == exact_divisor and exact_divisor.degree() == 121
    check_fibers(exact_divisor)
    assert conic_profile(SupplementalPole) == reference
    print("Both remaining cubic fibers match all 308 exact profiles.", flush=True)


class RecordedTangent(ExactRamifiedPole):
    def unit(self):
        if not super().unit():
            return False
        for coefficient in (self.coeffs[0], self.den):
            Recorded.record(tuple(int(c) % 3 for c in coefficient))
        return True


def ramified():
    reset_divisor()
    coordinates, exact = compute(RecordedTangent, 17, exact=True, last=12)
    assert len(exact) == 292
    exact_divisor = Recorded.divisor
    print("Exact tangent pole lattices through twelve.", flush=True)
    reset_divisor()
    residue_coordinates, residue = compute(RamifiedPole, 17)
    assert coordinates == residue_coordinates and exact == residue[:len(exact)]
    assert Recorded.divisor == exact_divisor and exact_divisor.degree() == 897
    verify_tangents(exact_divisor)
    assert compute(Supplemental, 17) == (coordinates, residue)
    print("All noncollision tangent fibers lift, with unit elimination at thirteen.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chart", nargs="?", choices=("nonbranch", "ramified"))
    chart = parser.parse_args().chart
    if chart != "ramified":
        nonbranch()
    if chart != "nonbranch":
        ramified()
