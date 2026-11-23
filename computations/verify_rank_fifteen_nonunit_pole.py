"""Thirty-one mixed-characteristic digits on the nonunit pole ray nu(y)=5 nu(x)>0."""

from flint import ctx, fmpq_poly, fmpz_poly, nmod_poly
import sympy as sp

from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_pole_fibers import RecordingExact, RecordingPole, cuts, tables
from verify_rank_fifteen_pole_slope_annulus import initial_tables, relations
from verify_rank_fifteen_unit_pole_depth import conic_profile
from verify_rank_fourteen_all_ordinary import gaussian_jets


def first_ray():
    v, m = sp.symbols("v m")
    l = 4*(1-v)**2/(1+v)**4
    x = v*(4+m*m*v)*l
    y = m**4*v**5*l**3
    slope = m**4/((4+m*m*v)**5*l*l)
    assert sp.cancel(y-x**5*slope) == 0
    assert sp.diff(x, v).subs(v, 0) == 16
    assert slope.subs(v, 0) == m**4/16384
    assert sp.diff(slope, m).subs(v, 0) == m**3/4096

    orders = []

    class Nonunit(RecordingExact):
        record_pivot_orders = staticmethod(orders.append)
        modulus = fmpz_poly([-1, 1])

    def local(field, half, s, scale):
        return relations(field, half, scale*s, s)

    Recorded.divisor = nmod_poly([1], 3)
    Recorded.record.cache_clear()
    products, coordinates = tables(Nonunit, 31, exact=True, build=local,
                         uncut_losses=(0, 1, 1, 4, 4), fifth_loss=6)
    cuts(Nonunit, products, losses=(4, 3, 1, 3, 1, 0, 1, 0))
    assert Nonunit.precision == 1 and Recorded.divisor == nmod_poly([0, 1], 3)
    assert len(orders) == 308
    return coordinates, orders


def outer_diagonal():
    precision = 47
    ctx.cap = precision+2
    jet = gaussian_jets(-1, precision, residue_prime=3)
    orders = []
    jet.record_pivot_orders = staticmethod(orders.append)

    def local(field, half, s, unused):
        return relations(field, half, s*s, s)

    products, coordinates = initial_tables(
        jet, precision, jet(fmpq_poly([1])/2), jet(fmpq_poly([0, 1])), jet.constant(1),
        build=local, uncut_losses=(0, 2, 2, 6, 6), center_losses=(0, 2, 2, 2), fifth_loss=8,
    )
    cuts(jet, products, losses=(6, 4, 2, 4, 2, 0, 2, 0))
    assert jet.precision == 1 and len(orders) == 308
    return coordinates, orders


def verify():
    ray_coordinates, ray = first_ray()
    boundary_coordinates, boundary = conic_profile(RecordingPole)
    diagonal_coordinates, diagonal = outer_diagonal()
    assert ray_coordinates == boundary_coordinates == diagonal_coordinates
    for a, b, c in zip(ray, boundary, diagonal, strict=True):
        assert tuple(i+j for i, j in zip(a, b, strict=True)) == tuple(c)
    print("Verified the nonunit pole sector 0 < 5 nu(x) <= nu(y).")


if __name__ == "__main__":
    verify()
