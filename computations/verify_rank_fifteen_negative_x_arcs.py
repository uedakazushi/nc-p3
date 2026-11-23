"""Two mixed-characteristic neighborhoods in the open negative-x cone."""

from flint import ctx, fmpq_poly

from verify_rank_fifteen_negative_x_boundary import diagonal, negative_x_relations
from verify_rank_fifteen_pole_fibers import RecordingPole, cuts
from verify_rank_fifteen_pole_slope_annulus import initial_tables
from verify_rank_fifteen_unit_pole_depth import conic_profile
from verify_rank_fourteen_all_ordinary import gaussian_jets


def curve(a, b, precision, uncut, centers, fifth, overlaps):
    ctx.cap = precision+2
    field = gaussian_jets(-1, precision, residue_prime=3)
    orders = []
    field.record_pivot_orders = staticmethod(orders.append)
    w = field(fmpq_poly([0]*a+[1]))
    g = field(fmpq_poly([0]*b+[1]))
    products, frame = initial_tables(
        field, precision, field(fmpq_poly([1])/2), w, g, build=negative_x_relations,
        uncut_losses=uncut, center_losses=centers, fifth_loss=fifth,
    )
    cuts(field, products, losses=overlaps)
    assert field.precision == 1 and len(orders) == 308
    return frame, orders


def shallow_y():
    return curve(1, 2, 37, (0, 2, 2, 4, 4), (0, 2, 4, 4), 4, (6, 2, 2, 2, 2, 0, 2, 0))


def shallow_g():
    return curve(2, 1, 21, (0, 1, 1, 3, 2), (0, 1, 2, 2), 2, (4, 1, 1, 1, 1, 0, 1, 0))


def verify():
    frame, middle = diagonal()
    edge_frame, edge = conic_profile(RecordingPole)
    y_frame, unit_y = shallow_y()
    g_frame, other = shallow_g()
    assert frame == edge_frame == y_frame == g_frame
    predicted = [tuple(x+y for x, y in zip(a, b, strict=True))
                 for a, b in zip(middle, edge, strict=True)]
    assert sum(a == b for a, b in zip(predicted, unit_y, strict=True)) == 254
    assert sum(a != b for a, b in zip(predicted[48:64], unit_y[48:64], strict=True)) == 4
    assert [j for j, (a, b) in enumerate(zip(middle, other, strict=True)) if a != b] == [81]
    assert [(a, b) for a, b in zip(middle[81], other[81], strict=True) if a != b] == [(0, 1)]*3
    print("Verified the two negative-x valuation tubes and all 308 comparisons.")


if __name__ == "__main__":
    verify()
