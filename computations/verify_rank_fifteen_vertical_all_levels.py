"""Descend the four central cuts from the vertical elliptic cover in characteristic three."""

from functools import cache

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx, nmod_poly

from verify_rank_fifteen_unit_residues import element, verify_unit_point

from rank_fourteen_function_field import CurveElement
from verify_rank_fifteen_cubic_zero_disk import tensor_relations
from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_parameter_orbits import pattern
from verify_rank_fifteen_pencil_trace_disks import local_central_sections
from verify_rank_fifteen_pole_fibers import transform
from verify_rank_fourteen_all_ordinary import DIMENSIONS, saturated_projection
from verify_rank_fourteen_central_cut import quotient_layers

ZERO = nmod_poly([], 3)
ONE = nmod_poly([1], 3)
T = nmod_poly([0, 1], 3)


class EllipticUnit(CurveElement):
    """Coordinates in the basis 1, i, W, iW of F3(i,t,W) over F3(t)."""

    radicand = T-T**3
    precision = 1

    def __init__(self, coordinates, denominator=ONE):
        assert denominator
        if not any(coordinates):
            self.coordinates, self.denominator = (ZERO,)*4, ONE
            return
        if denominator != ONE:
            common = denominator
            for coefficient in coordinates:
                if common.degree() == 0:
                    break
                common = common.gcd(coefficient)
            if common.degree() > 0:
                coordinates = tuple(c//common for c in coordinates)
                denominator = denominator//common
            if denominator.leading_coefficient() != 1:
                coordinates = tuple(-c for c in coordinates)
                denominator = -denominator
        self.coordinates, self.denominator = coordinates, denominator

    @classmethod
    def lift(cls, value):
        if isinstance(value, cls):
            return value
        coefficient = value if isinstance(value, nmod_poly) else nmod_poly([int(value)%3], 3)
        return cls((coefficient, ZERO, ZERO, ZERO))

    constant = lift

    def unit(self):
        return bool(self)

    def inverse(self):
        Recorded.record(tuple(map(int, self.norm_numerator())))
        return super().inverse()


def windows():
    t = EllipticUnit.lift(T)
    i = EllipticUnit((ZERO, ONE, ZERO, ZERO))
    w = EllipticUnit((ZERO, ZERO, ONE, ZERO))
    a, b = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    v, square_root = 2*a*b, w/(1+t*t)
    assert square_root**2 == v
    k = 1/v
    h, n = k-v, (a*a-b*b)*a/v
    p, u, q, r = (v+1)/square_root, i*(a+b), (v-1)/square_root, a-b
    assert (p*p, u*u, q*q, r*r) == (2*(k+1+h), 2*(k+1-h), 2*(k-1+h), 2*(k-1-h))
    assert n*n+(h/n)**2 == k*k-1
    one = EllipticUnit.constant(1)
    gamma = pattern(one, (k-1)/2, n/2, h/(2*n), (u-p)/4, (u+p)/4, (r-q)/4, (r+q)/4)
    return tensor_relations(gamma, EllipticUnit.constant(0))


def sparse_pivot(rows):
    def score(position):
        j, k = position
        entry = rows[j][k]
        degree = max(c.degree() for c in (*entry.coordinates, entry.denominator))
        return len(rows[j])//3, degree, len(rows[j]), k
    return min(((j, k) for j, row in enumerate(rows) for k in row), key=score)


def pivot_divisor(build_windows=windows):
    Recorded.divisor = ONE
    Recorded.record.cache_clear()
    relations = build_windows()

    def project(rows, width, target, field):
        return saturated_projection(rows, width, target, field, select_pivot=sparse_pivot)

    centers = local_central_sections(relations, EllipticUnit, project=project,
                                     invariant_pairs=((0, 1, 0), (0, 0, 1)))

    def cut(n, h, rows, width):
        free, images, loss = project(rows, width, DIMENSIONS[n][h], EllipticUnit)
        assert loss == 0
        return free, images

    for start in range(4):
        for degree in quotient_layers(relations, centers, start, 13, EllipticUnit.constant(1), cut):
            print("Start", start, "length", degree, "divisor degree", Recorded.divisor.degree(), flush=True)
    return Recorded.divisor


def homogeneous_pullback(f, a, b):
    """Evaluate the homogeneous substitution without quadratic Horner growth."""
    @cache
    def powers(n):
        return a**n, b**n

    def visit(part, degree):
        if degree < 48:
            return transform(part, a, b)*b**(degree-part.degree()) if part else 0*a
        split = (degree+1)//2
        return (powers(degree-split+1)[1]*visit(part.truncate(split), split-1)
                + powers(split)[0]*visit(part.right_shift(split), degree-split))
    return visit(f, f.degree())


def residue_polynomial(value, modulus):
    orbit = [value]
    conjugate = value.pow_mod(3, modulus)
    while conjugate != value:
        orbit.append(conjugate)
        conjugate = conjugate.pow_mod(3, modulus)
    coefficients = [ONE]
    for root in orbit:
        result = [ZERO]*(len(coefficients)+1)
        for index, coefficient in enumerate(coefficients):
            result[index] = (result[index]-coefficient*root) % modulus
            result[index+1] = (result[index+1]+coefficient) % modulus
        coefficients = result
    assert all(f.degree() <= 0 for f in coefficients)
    return tuple(int(f[0]) for f in coefficients)


def descend(divisor):
    a, b = T*(1-T**2), (1+T**2)**2
    num, den = (a*b*b-a**3)**2, b**6
    assert den.degree() == 24 and num.gcd(den).degree() == 0
    # Away from y=0 and infinity the map is etale.
    assert not (num*den)**2 % (num.derivative()*den-num*den.derivative())
    core = divisor // divisor.gcd(num*den)
    assert core.gcd(core.derivative()).degree() == 0
    for c, d in ((-T, ONE), (ONE, T), (1+T, 1-T)):
        assert transform(num, c, d)*d**(24-num.degree())*den == transform(den, c, d)*num
    while core.degree() > 0:
        degree = core.degree()
        reflected = nmod_poly([int(c)*(-1 if j % 2 else 1) for j, c in enumerate(core)], 3)
        core = core.gcd(reflected)
        core = core.gcd(core.reverse())
        core = core.gcd(homogeneous_pullback(core, 1+T, 1-T))
        if core.degree() == degree:
            break

    groups = {}
    for modulus, exponent in core.factor()[1]:
        assert exponent == 1
        common, inverse, _ = den.xgcd(modulus)
        assert common == ONE
        polynomial = residue_polynomial(num*inverse % modulus, modulus)
        groups[polynomial] = groups.get(polynomial, ONE)*modulus

    exceptions = []
    for polynomial, preimage in sorted(groups.items()):
        fiber = transform(nmod_poly(list(polynomial), 3), num, den)
        fiber *= pow(int(fiber.leading_coefficient()), -1, 3)
        assert fiber.degree() == 24*(len(polynomial)-1) and not fiber % preimage
        if preimage == fiber:
            assert not divisor % fiber
            exceptions.append(polynomial)
    return core.degree(), tuple(exceptions)


def check_cover():
    divisor = pivot_divisor()
    assert divisor.degree() == 204451
    core, exceptional = descend(divisor)
    print("Stable divisor degree", core, "exceptional fibers", exceptional, flush=True)
    assert core == 784
    assert exceptional == ((1, 1), (1, 2, 1, 1), (1, 2, 2, 1, 0, 0, 1), (2, 1))


def check_sextic():
    ring = fmpz_mod_poly_ctx(3)
    modulus = ring([2, 2, 1, 0, 2, 0, 1])
    g = [1, 2, 2, 1, 0, 0, 1]
    assert modulus.is_irreducible() and ring(g).is_irreducible()
    field = fq_default_ctx(modulus=modulus, var="a")
    y = element(field, 124)
    assert fq_default_poly_ctx(field)(g)(y).is_zero()
    verify_unit_point(field, field.one(), y, (299, 378, 183, 137, 9, 500, 467), ((0, 1, 0), (0, 0, 1)))
    print("Verified all four starts on the exceptional sextic fiber.", flush=True)


if __name__ == "__main__":
    check_sextic()
    check_cover()
