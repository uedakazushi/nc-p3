"""The remaining nonunit conic sector and its generic s=0 boundary."""

from itertools import product

from flint import ctx, fmpq_poly, fmpz_poly, fq_default_ctx

from verify_rank_fifteen_nonunit_pole import first_ray
from verify_rank_fifteen_pole_fibers import cuts, tables
from verify_rank_fifteen_pole_slope_annulus import (
    IntegralPoleJet, initial_tables, pattern, relations, tensor_relations,
)
from verify_rank_fourteen_all_ordinary import gaussian_jets


def boundary_relations(field, half, s, q):
    one = field.constant(1)
    conic_inverse = (one+s*s).inverse()
    a, b = (one-s*s)*conic_inverse, 2*s*conic_inverse
    slope = a*b

    v = 2*q*(one+q*q).inverse()
    scale = (one+q*q)*(one-q*q).inverse()
    r, k = 2*v*(one-v*v).inverse(), (one+v*v)*(one-v*v).inverse()

    def root(square, initial):
        value = initial
        for _ in range(field.precision.bit_length()):
            value = (value+square*value.inverse())*half
        assert not (value*value-square)
        return value

    p, u = (root(2*(k+one)+sign*2*slope*r*r, 2*scale) for sign in (1, -1))
    small = scale.inverse()
    square = 2*(k+one).inverse()
    v, w = (r*root(square+sign*2*slope, sign*small) for sign in (1, -1))
    gamma = pattern(one, (k-one)*half, r*a*half, r*b*half, (u-p)*half*half,
                    (u+p)*half*half, (w-v)*half*half, (w+v)*half*half)
    return tensor_relations(gamma, one*0)


def boundary():
    orders = []

    class Field(IntegralPoleJet):
        modulus = fmpz_poly([1, 0, 2, 1])
        record_pivot_orders = staticmethod(orders.append)

    products, coordinates = tables(Field, 15, exact=True, build=boundary_relations,
                                  uncut_losses=(0, 0, 0, 2, 2),
                                  center_losses=(0, 0, 0, 0), fifth_loss=4)
    cuts(Field, products, losses=(2, 2, 0, 2, 0, 0, 0, 0))
    assert Field.precision == 1 and len(orders) == 308
    return coordinates, orders


def diagonal():
    ctx.cap = 47
    jet = gaussian_jets(-1, 45, residue_prime=3)
    orders = []
    jet.record_pivot_orders = staticmethod(orders.append)

    def build(field, half, t, unused):
        return relations(field, half, t, t*t)

    products, coordinates = initial_tables(
        jet, 45, jet(fmpq_poly([1])/2), jet(fmpq_poly([0, 1])), jet.constant(1),
        build=build, uncut_losses=(0, 1, 1, 6, 6), fifth_loss=10,
    )
    cuts(jet, products, losses=(6, 5, 1, 5, 1, 0, 1, 0))
    assert jet.precision == 1 and len(orders) == 308
    return coordinates, orders


def binary_rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            column = value.bit_length()-1
            if column not in pivots:
                pivots[column] = value
                break
            value ^= pivots[column]
    return len(pivots)


def unit_cover():
    field = fq_default_ctx(3, 2, "a")
    one = field.one()
    i = (-one).sqrt()
    points = [field([a, b]) for a, b in product(range(3), repeat=2)]

    def supports(sign, shift):
        masks = [1 << 9 for _ in range(6)]  # Every cubic has an odd pole at infinity.
        for j, old in enumerate(points):
            v = sign*old+shift
            radicands = []
            for epsilon in (1, -1):
                radicands.extend(((one-v*v)*(one+epsilon*2*i*v),
                                  v*(one+v)*(one+epsilon*i*(one-v)),
                                  v*(one-v)*(one+epsilon*i*(one+v))))
            for k, radicand in enumerate(radicands):
                if radicand.is_zero():
                    masks[k] |= 1 << j
        assert all(mask.bit_count() == 4 for mask in masks)
        return masks

    basis = supports(1, 0)
    assert binary_rank(basis) == 6
    for j, v in enumerate(points):
        assert any(mask & (1 << j) for mask in basis)
        y = -(v-v**3)**2
        assert y == (field.zero() if v**3 == v else one)
    for sign, shift in product((1, -1), range(3)):
        assert binary_rank(basis+supports(sign, shift)) == 6
    print("Verified the genus-97 signed cover of degree 384 on xbar=-1.", flush=True)


def verify():
    unit_cover()
    edge_coordinates, edge = boundary()
    ray_coordinates, ray = first_ray()
    middle_coordinates, middle = diagonal()
    assert edge_coordinates == ray_coordinates == middle_coordinates
    for a, b, c in zip(edge, ray, middle, strict=True):
        assert tuple(j+k for j, k in zip(a, b, strict=True)) == tuple(c)
    print("Verified every point with positive x valuation and 3 nu(x) < nu(y) < 5 nu(x).")


if __name__ == "__main__":
    verify()
