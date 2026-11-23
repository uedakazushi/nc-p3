"""Match five collision boundary/diagonal charts in common central coordinates."""

from functools import cache

from flint import ctx, fmpq_poly, fmpz_poly, nmod_poly
import sympy as sp

from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_pole_fibers import RecordingExact, cuts, tables
from verify_rank_fifteen_pole_slope_annulus import initial_tables
from verify_rank_fifteen_ramified_pole import relations
from verify_rank_fourteen_all_ordinary import gaussian_jets, saturated_projection


def fixed_frame(target):
    def frame(h, basis, image, jet):
        words_by_degree = basis.copy()
        wanted = set()

        class BasisCheck:
            precision = jet.precision
            constant = staticmethod(jet.constant)

        for character, start, words, indices in target:
            if h != character:
                continue
            if list(words) != basis[start, 4, h]:
                rows = [image(start, word).copy() for word in words]
                assert saturated_projection(rows, len(words), 0, BasisCheck)[2] == 0
                words_by_degree[start, 4, h] = list(words)
            offset = sum(len(basis[p, 4, h]) for p in range(start))
            wanted.update(offset+i for i in indices)

        def select(rows):
            return next(((i, k) for i in sorted(range(len(rows)), key=lambda j: len(rows[j]))
                         for k in sorted(rows[i]) if k not in wanted and rows[i][k].unit()), None)
        return words_by_degree, select
    return frame


PROFILES = {
    (1, 1): ((0, 2, 2, 3, 2), (0, 0, 0, 0), 0, 0, (2, 0, 1, 0, 1, 0, 1, 0)),
    (1, 2): ((0, 2, 2, 3, 2), (1, 0, 0, 0), 0, 0, (2, 0, 1, 0, 1, 0, 1, 0)),
    (2, 1): ((0, 4, 4, 5, 4), (0, 1, 1, 1), 0, 1, (3, 0, 2, 1, 2, 0, 2, 0)),
    (3, 1): ((0, 6, 6, 7, 6), (0, 2, 2, 2), 0, 2, (4, 0, 3, 2, 4, 0, 3, 0)),
    (3, 2): ((0, 6, 6, 8, 6), (0, 1, 1, 1), 0, 1, (5, 0, 3, 1, 3, 0, 3, 0)),
    (2, 3): ((0, 4, 4, 6, 4), (1, 0, 0, 0), 0, 0, (4, 0, 2, 0, 2, 0, 2, 0)),
    (1, 0): ((0, 2, 2, 2, 2), (0, 1, 1, 1), 0, 1, (1, 0, 1, 1, 2, 0, 1, 0)),
}


def profile(weights, *, exact=False, target=None, modulus=(-1, 1)):
    a, b = weights
    uncut, center, fourth, fifth, overlaps = PROFILES[weights]
    precision = sum(uncut)+max(center)+fifth+sum(overlaps)+1
    options = dict(uncut_losses=uncut, center_losses=center, fourth_loss=fourth, fifth_loss=fifth,
                   center_frame=fixed_frame(target) if target else None)
    orders = []

    def build(field, half, t, m):
        def power(n):
            value = field.constant(1)
            for _ in range(n):
                value = value*t
            return value
        if b == 0:
            return relations(field, half, t, m)
        h, delta = power(a), power(b)
        if a > b:
            h = h*m
        else:
            delta = delta*m
        return relations(field, half, h, field.constant(1)+delta)

    if exact:
        class Field(RecordingExact):
            record_pivot_orders = staticmethod(orders.append)
        Field.modulus = fmpz_poly(list(modulus))
        Recorded.divisor = nmod_poly([1], 3)
        Recorded.record.cache_clear()
        products, coordinates = tables(Field, precision, exact=True, build=build, **options)
    else:
        ctx.cap = precision+2
        Field = gaussian_jets(-1, precision, residue_prime=3)
        Field.record_pivot_orders = staticmethod(orders.append)
        products, coordinates = initial_tables(Field, precision, Field(fmpq_poly([1])/2),
                                              Field(fmpq_poly([0, 1])), Field.constant(1),
                                              build=build, **options)
    cuts(Field, products, losses=overlaps)
    assert Field.precision == 1 and len(orders)==308
    if target:
        assert coordinates==target
    return coordinates, orders


def compare(first, second, middle):
    c1, left = first
    c2, right = second
    cm, diagonal = middle
    assert c1 == c2 == cm
    for a, b, c in zip(left, right, diagonal, strict=True):
        assert tuple(i+j for i, j in zip(a, b, strict=True)) == tuple(c)


@cache
def unit_ray():
    primary = profile((2, 1), exact=True)
    t = nmod_poly([0, 1], 3)
    assert Recorded.divisor == t*(t+1)
    profile((2, 1), exact=True, modulus=(1, 1))
    assert Recorded.divisor == t*(t-1)
    return primary


def shallow_cone():
    middle = profile((3, 1))
    boundary = profile((1, 0), exact=True, target=middle[0], modulus=(2, 1, 1))
    compare(boundary, unit_ray(), middle)
    print("Matched all 308 lists on the shallow collision cone.", flush=True)


def other_cones():
    compare(unit_ray(), profile((1, 1), exact=True), profile((3, 2)))
    middle = profile((2, 3))
    aligned = profile((1, 1), exact=True, target=middle[0])
    deepest_ray = profile((1, 2), exact=True, modulus=(0, 1))
    compare(aligned, deepest_ray, middle)
    print("Matched all 616 lists on both remaining collision cones.", flush=True)


def vertical_levels():
    a, t = sp.symbols("a t", nonzero=True)
    f = (1-a)**2-4*a**3/27
    y = a**3/27-(t-f/t)**2/16
    p, q = 3*sp.sqrt(3)*(t-f/t)/4, (t+f/t)/2
    assert sp.cancel(p*p-a**3+27*y) == 0 and sp.cancel(q*q-(1-a)**2+4*y) == 0
    assert sp.cancel(f+(a-3)**2*(4*a-3)/27) == 0
    assert sp.cancel(q+2*p/(3*sp.sqrt(3))-t) == 0
    exceptional = (3+sp.sqrt(2))/2
    assert sp.expand(4*exceptional**2-12*exceptional+7) == 0
    assert (4-12+7) % 3 != 0


if __name__ == "__main__":
    vertical_levels()
    shallow_cone()
    other_cones()
