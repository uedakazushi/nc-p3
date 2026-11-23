"""Fold the two conic poles to one square coordinate and certify the limiting jet."""

from functools import partial

import sympy as sp
from flint import ctx, fmpq_poly

from verify_rank_fifteen_negative_x_arcs import shallow_g
from verify_rank_fifteen_negative_x_boundary import negative_x_relations
from verify_rank_fifteen_pole_fibers import cuts
from verify_rank_fifteen_pole_slope_annulus import initial_tables, pattern, tensor_relations
from verify_rank_fourteen_all_ordinary import gaussian_jets


def folded_relations(field, half, z, u, *, imaginary_unit):
    one = field.constant(1)
    i = imaginary_unit

    def root(square, initial):
        value = initial
        for _ in range(field.precision.bit_length()):
            value = (value+square*value.inverse())*half
        assert not (value*value-square)
        return value

    v = -i*root(one-z, one)
    k = root(one+z*u*u, one)
    p, other_p = (root(2*(k+one+sign*v*u*u), 2*one) for sign in (1, -1))
    q, other_q = (u*root(2*z*(k+one).inverse()+sign*2*v, one-sign*i)
                  for sign in (1, -1))
    gamma = pattern(one, (k-one)*half, v*u*half, u*half,
                    (other_p-p)*half*half, (other_p+p)*half*half,
                    (other_q-q)*half*half, (other_q+q)*half*half)
    return tensor_relations(gamma, field.constant(0))


def verify():
    w, g = sp.symbols("w g")
    s = sp.I+w
    u, h = 2*s*g/(2*sp.I+w), (1+s*s)/(2*s)
    v = (1-s*s)/(2*s)
    assert sp.cancel(h*u-g*w) == 0 and sp.cancel(v*v-h*h+1) == 0
    assert sp.cancel(u*v-g*(1-s*s)/(2*sp.I+w)) == 0
    z, r, slope = h*h, w*g, v/(h*h)
    assert sp.cancel((4+r*r)*slope**2-(4+z*u*u)*(z-1)/z**2) == 0
    assert sp.cancel(r**4*slope**6-u**4*(z-1)**3/z**4) == 0

    ctx.cap = 23
    field = gaussian_jets(-1, 21, residue_prime=3)
    t = field(fmpq_poly([0, 1]))
    i = field(fmpq_poly([]), fmpq_poly([1]))
    half = field(fmpq_poly([1])/2)
    zero = field.constant(0)
    build = partial(folded_relations, imaginary_unit=i)
    source, target = negative_x_relations(field, half, zero, t), build(field, half, zero, t)
    for side, other_side in zip(source, target, strict=True):
        for (h1, row), (h2, other_row) in zip(side, other_side, strict=True):
            assert h1 == h2 and row.keys() == other_row.keys()
            assert all(not (v-other_row[word]) for word, v in row.items())

    orders = []
    field.record_pivot_orders = staticmethod(orders.append)
    products, frame = initial_tables(field, 21, half, zero, t, build=build,
                        uncut_losses=(0, 1, 1, 3, 2), center_losses=(0, 1, 2, 2), fifth_loss=2)
    cuts(field, products, losses=(4, 1, 1, 1, 1, 0, 1, 0))
    assert field.precision == 1 and len(orders) == 308
    assert shallow_g() == (frame, orders)
    print("Verified the folded pole coordinate and the full sector 10 nu(y) <= 19 nu(x) < 0.")


if __name__ == "__main__":
    verify()
