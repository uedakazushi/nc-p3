"""Simultaneous two-parameter Smith elimination on three outer pole cones."""

from flint import ctx, fmpq_poly, fmpz_poly, nmod_poly

from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_outer_pole import outer_relations
from verify_rank_fifteen_pole_fibers import RecordingExact, RecordingPole, cuts, tables
from verify_rank_fifteen_pole_slope_annulus import initial_tables
from verify_rank_fifteen_ramified_pole import RamifiedPole, compute
from verify_rank_fifteen_unit_pole_depth import conic_profile
from verify_rank_fourteen_all_ordinary import gaussian_jets

# (ord w, ord u): successive precision losses on uncut paths, the center,
# degree five of the cut, and the two-sided overlap matrices at lengths 6--13.
PROFILES = {
    (1, 1): ((0, 3, 3, 4, 4), (0, 2, 2, 2), 3, (3, 1, 2, 1, 3, 0, 2, 0)),
    (1, 2): ((0, 4, 4, 6, 6), (0, 3, 3, 3), 5, (5, 2, 3, 2, 4, 0, 3, 0)),
    (2, 1): ((0, 5, 5, 6, 6), (0, 3, 3, 3), 4, (4, 1, 3, 2, 5, 0, 3, 0)),
    (2, 3): ((0, 7, 7, 10, 10), (0, 5, 5, 5), 8, (8, 3, 5, 3, 7, 0, 5, 0)),
    (1, 3): ((0, 5, 5, 8, 8), (0, 4, 4, 4), 7, (7, 3, 4, 3, 5, 0, 4, 0)),
}


def diagonal(weights):
    uncut, center, fifth, overlaps = PROFILES[weights]
    precision = sum(uncut)+max(center)+fifth+sum(overlaps)+1
    ctx.cap = precision+2
    jet = gaussian_jets(-1, precision, residue_prime=3)
    orders = []
    jet.record_pivot_orders = staticmethod(orders.append)

    def build(field, half, t, unused):
        def power(n):
            value = field.constant(1)
            for _ in range(n):
                value = value*t
            return value
        return outer_relations(field, half, power(weights[0]), power(weights[1]))

    products, coordinates = initial_tables(
        jet, precision, jet(fmpq_poly([1])/2), jet(fmpq_poly([0, 1])), jet.constant(1),
        build=build, uncut_losses=uncut, center_losses=center, fifth_loss=fifth,
    )
    cuts(jet, products, losses=overlaps)
    assert jet.precision == 1 and len(orders) == 308
    return coordinates, orders


def check_second_ray():
    parameter = nmod_poly([0, 1], 3)
    Recorded.divisor = 1+0*parameter
    Recorded.record.cache_clear()
    orders = []

    class OuterRay(RecordingExact):
        modulus = fmpz_poly([-1, 1])
        record_pivot_orders = staticmethod(orders.append)

    def build(jet, half, w, m):
        return outer_relations(jet, half, w, m*w*w)

    uncut, center, fifth, overlaps = PROFILES[1, 2]
    products, coordinates = tables(
        OuterRay, 48, exact=True, build=build, uncut_losses=uncut,
        center_losses=center, fifth_loss=fifth,
    )
    cuts(OuterRay, products, losses=overlaps)
    assert OuterRay.precision == 1 and Recorded.divisor == parameter and len(orders) == 308
    print("All unit residues lift on u=m w^2.", flush=True)
    return coordinates, orders


def check_cones():
    charts = {(1, 2): check_second_ray()}
    charts[0, 1] = conic_profile(RecordingPole)
    charts[1, 0] = compute(RamifiedPole, 17)
    for weights in PROFILES:
        if weights not in charts:
            charts[weights] = diagonal(weights)
    for left, right in (((1, 0), (1, 1)), ((1, 1), (1, 2)), ((1, 2), (0, 1))):
        middle = tuple(a+b for a, b in zip(left, right, strict=True))
        edge_coordinates, first = charts[left]
        other_coordinates, second = charts[right]
        diagonal_coordinates, check = charts[middle]
        assert edge_coordinates == other_coordinates == diagonal_coordinates
        for a, b, c in zip(first, second, check, strict=True):
            assert tuple(j+k for j, k in zip(a, b, strict=True)) == tuple(c)
        print("Matched all 308 Smith lists on the cone", left, right, flush=True)


if __name__ == "__main__":
    check_cones()
