"""Three rays and their open cones in the square coordinate of the negative-x chart."""

from functools import partial

from flint import fmpz_poly, fq_default_ctx, fq_default_poly_ctx, nmod_poly

from rank_fifteen_gaussian_parameter import gaussian_parameter_jets
from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_negative_x_fold import folded_relations
from verify_rank_fifteen_pole_fibers import RecordingPole, cuts
from verify_rank_fifteen_pole_slope_annulus import initial_tables
from verify_rank_fifteen_unit_pole_depth import conic_profile


def ray_checks(d):
    if d == 1:
        return 19, (0, 1, 1, 2, 2), (0, 1, 2, 2), 2, (3, 1, 1, 1, 1, 0, 1, 0)
    return 21, (0, 1, 1, 3, 2), (0, 1, 2, 2), 2, (4, 1, 1, 1, 1, 0, 1, 0)


def tables(field, z, u, i, half, checks):
    digits, uncut, center, fifth, overlaps = checks
    orders = []
    field.record_pivot_orders = staticmethod(orders.append)
    products, frame = initial_tables(field, digits, half, z, u,
                        build=partial(folded_relations, imaginary_unit=i),
                        uncut_losses=uncut, center_losses=center, fifth_loss=fifth)
    cuts(field, products, losses=overlaps)
    assert field.precision == 1 and len(orders) == 308
    return frame, orders


def generic(d, modulus):
    checks = ray_checks(d)
    field = gaussian_parameter_jets(modulus)
    field.precision = checks[0]
    base = field.base
    zero, one, m = fmpz_poly([]), fmpz_poly([1]), fmpz_poly([0, 1])
    i = field(base([]), base([one]))
    z, u = field(base([zero]*d+[m])), field(base([zero, one]))
    half = field(base([one], 2*one))
    Recorded.record.cache_clear()
    Recorded.divisor = nmod_poly([1], 3)
    profile = tables(field, z, u, i, half, checks)
    return profile, Recorded.divisor


def residue_jets(k):
    ring = fq_default_poly_ctx(k)

    class Jet:
        __slots__ = ("value",)

        def __init__(self, value):
            self.value = value

        @classmethod
        def constant(cls, value):
            return cls(ring(value))

        def __bool__(self):
            return bool(self.value)

        def __neg__(self):
            return Jet(-self.value)

        def __add__(self, other):
            right = other.value if isinstance(other, Jet) else ring(other)
            return Jet(self.value.add_trunc(right, self.precision))

        __radd__ = __add__

        def __sub__(self, other):
            return self + -other

        def __rsub__(self, other):
            return -self + other

        def __mul__(self, other):
            right = other.value if isinstance(other, Jet) else ring(other)
            return Jet(self.value.mul_low(right, self.precision))

        __rmul__ = __mul__

        def unit(self):
            return not self.value[0].is_zero()

        def inverse(self):
            assert self.unit()
            return Jet(self.value.inverse_series_trunc(self.precision))

        def __truediv__(self, other):
            return self * other.inverse()

        def divide_by_parameter(self):
            assert not self.unit()
            return Jet(self.value.right_shift(1))

    return ring, Jet


def specialize(k, value, z_power, u_power, checks):
    ring, field = residue_jets(k)
    t = ring.gen()
    field.precision = checks[0]
    z, u = field(t**z_power*value), field(t**u_power)
    i = field(ring([(-k.one()).sqrt()]))
    return tables(field, z, u, i, field.constant(2), checks)


def compare_cone(k, z_power, u_power, bounds, first, second):
    frame, current = specialize(k, k.one(), z_power, u_power, bounds)
    left_frame, left = first
    right_frame, right = second
    assert frame == left_frame == right_frame
    for a, b, c in zip(left, right, current, strict=True):
        assert tuple(i+j for i, j in zip(a, b, strict=True)) == c


def check_residues(d, divisor, generic_profile):
    # The two bounding rays d=1,2 require unit m. The last chart allows m=0 as well.
    for factor, exponent in divisor.factor()[1]:
        assert exponent == 1
        if d < 3 and factor == nmod_poly([0, 1], 3):
            continue
        degree = factor.degree()
        extension = degree if degree % 2 == 0 else 2*degree
        k = fq_default_ctx(3, extension, "a")
        roots = fq_default_poly_ctx(k)(list(map(int, factor))).roots()
        assert len(roots) == degree
        for residue, multiplicity in roots:
            assert multiplicity == 1
            assert specialize(k, residue, d, 1, ray_checks(d)) == generic_profile


def verify():
    k = fq_default_ctx(3, 2, "a")
    curves = []
    for d, modulus in ((1, (1, 2, 0, 1)), (2, (2, 1)), (3, (1, 2, 0, 1))):
        profile, divisor = generic(d, modulus)
        m = nmod_poly([0, 1], 3)
        expected = {1: m*(m**4+1), 2: m*(m+1), 3: nmod_poly([1], 3)}
        assert divisor == expected[d]
        print("Ray", d, "pivot divisor", divisor.factor(), flush=True)
        check_residues(d, divisor, profile)
        curves.append(profile)

    edge = conic_profile(RecordingPole)
    overlaps = (2, 2, 2, 2, 0, 2, 0)
    compare_cone(k, 1, 2, (35, (0, 2, 2, 4, 4), (0, 2, 3, 3), 4, (5, *overlaps)),
                 edge, curves[0])
    compare_cone(k, 3, 2, (39, (0, 2, 2, 5, 4), (0, 2, 4, 4), 4, (7, *overlaps)),
                 curves[0], curves[1])
    compare_cone(k, 5, 2, (41, (0, 2, 2, 6, 4), (0, 2, 4, 4), 4, (8, *overlaps)),
                 curves[1], curves[2])
    print("Certified every ordinary rank-fifteen window with negative x valuation.")


if __name__ == "__main__":
    verify()
