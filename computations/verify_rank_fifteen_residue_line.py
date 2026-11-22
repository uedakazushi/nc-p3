"""Four central cuts on the entire x=0 geometric residue line in characteristic three."""

from concurrent.futures import ProcessPoolExecutor
from functools import cache

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx

from rank_fourteen_dual import fitting_support, polynomial
from verify_klein_center_lift import central_equations, echelon, field_paths
from verify_rank_fourteen_central_cut import central_sections
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile
from verify_rank_fifteen_parameter_orbits import amplitudes, dot, pattern
from verify_rank_fifteen_pencil_linear_residues import graded_profile
from verify_rank_fifteen_regular_orbits import tensor_windows

F9 = fq_default_ctx(modulus=fmpz_mod_poly_ctx(3)([1, 0, 1]), var="i")
R = fq_default_poly_ctx(F9)
ONE = R.one()


def coefficients(poly):
    return tuple(tuple(map(int, c.to_list())) for c in poly)


def from_coefficients(encoded):
    return R([F9(list(pair)) for pair in encoded])


def integral_gamma():
    """The normalized Pauli coefficients multiplied by the common unit 4s."""
    s = R.gen()
    p, q = s**2 + 1, s**2 - 1
    u, v = q * F9.gen(), p * F9.gen()
    return pattern(4 * s, -2 * s, p - u, p + u, u - p, u + p, v - q, v + q)


def integral_windows():
    gamma = integral_gamma()

    def component(a, b, c, d):
        h = a ^ c
        if h != b ^ d:
            return R.zero()
        return sum(
            (gamma[h + 4 * j] * (-1) ** (dot(j, c) + dot((0, 3, 2, 1)[j], d))
             for j in range(4)), R.zero()
        )

    windows = []
    for start in range(4):
        side = []
        for h in range(4):
            for k in range(4):
                row = {}
                for j in range(4):
                    word = (j, j ^ h, k, k ^ h)
                    value = component(*(word[-start:] + word[:-start]))
                    if not value.is_zero():
                        row[j, j ^ h] = value
                if row:
                    side.append((h, row))
        windows.append(side)
    return windows


class Fraction:
    """Arithmetic in F9(s), recording the numerators inverted in path and center charts."""

    __slots__ = ("n", "d")
    inverted = set()

    def __init__(self, numerator, denominator=ONE):
        if denominator.is_zero():
            raise ZeroDivisionError
        if not denominator.is_one():
            common = numerator.gcd(denominator) * denominator.leading_coefficient()
            numerator, denominator = numerator // common, denominator // common
        self.n, self.d = numerator, denominator

    @staticmethod
    @cache
    def constant(value):
        return Fraction(R([value]))

    def __bool__(self):
        return not self.n.is_zero()

    def __neg__(self):
        return Fraction(-self.n, self.d)

    def __add__(self, other):
        if not isinstance(other, Fraction):
            other = self.constant(other)
        if self.d == other.d:
            return Fraction(self.n + other.n, self.d)
        return Fraction(self.n * other.d + other.n * self.d, self.d * other.d)

    __radd__ = __add__

    def __mul__(self, other):
        if not isinstance(other, Fraction):
            other = self.constant(other)
        return Fraction(self.n * other.n, self.d * other.d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, Fraction):
            other = self.constant(other)
        if other.n.degree() > 0:
            self.inverted.add(coefficients(other.n.monic()))
        return Fraction(self.n * other.d, self.d * other.n)


def encode(poly):
    """Represent a polynomial as its real and imaginary F3[s] coordinates."""
    return tuple(polynomial([pair[k] for pair in coefficients(poly)]) for k in range(2))


def clear_row(row):
    denominator = ONE
    for value in row.values():
        denominator *= value.d // denominator.gcd(value.d)
    return {word: encode(value.n * (denominator // value.d)) for word, value in row.items()}


def path_and_center(windows):
    Fraction.inverted.clear()
    unit = Fraction.constant(1)
    rational = [[(h, {w: Fraction(c) for w, c in row.items()}) for h, row in side]
                for side in windows]
    basis, image = field_paths(rational, unit)
    for p in range(4):
        assert [sum(len(basis[p, n, h]) for h in range(4)) for n in range(6)] == [
            1, 4, 10, 20, 35, 56
        ]
    for h, expected in enumerate((3, 2, 0, 0)):
        rows, width = central_equations(basis, image, h)
        assert width - len(echelon(rows)) == expected
    sections = central_sections(basis, image, unit, ((1, 3, 5), (2, 7, 11)))
    cuts = [[(h, clear_row(row)) for h, row in side] for side in sections]
    bad = set()
    for pivot in Fraction.inverted:
        bad.update(coefficients(f) for f, _ in from_coefficients(pivot).factor()[1])
    return cuts, bad


def remaining_orbits(bad):
    """Only complete orbits of irreducibles can obstruct all choices above one y."""
    s, i = R.gen(), F9.gen()
    bad = set(bad)
    pending = []
    while bad:
        element = next(iter(bad))
        f = from_coefficients(element)
        if f == s or ((s**4 + 1) % f).is_zero():
            bad.remove(element)
            continue
        orbit = set()
        for u in (F9.one(), -F9.one(), i, -i):
            scaled = R([c * u**j for j, c in enumerate(f)])
            orbit.update(coefficients(q.monic()) for q in (scaled, scaled.reverse()))
        if orbit <= bad:
            pending.append(min(orbit))
        bad.difference_update(orbit)
    return sorted(pending)


def check_point(encoded):
    degree = len(encoded) - 1
    field = fq_default_ctx(3, 2 * degree, "a")
    i = field(-1).sqrt()
    poly = fq_default_poly_ctx(field)([field(a) + i*b for a, b in encoded])
    s = poly.roots()[0][0]
    n = (s + 1/s - i*(s - 1/s)) / 2
    gamma = amplitudes(s, i, n, field.one())
    for c, polynomial_value in zip(gamma, integral_gamma()):
        value = field.zero()
        for a, b in reversed(coefficients(R(polynomial_value))):
            value = value*s + a + i*b
        assert value == 4*s*c
    windows = tensor_windows(gamma, field)
    unit = Scalar(field.one())
    try:
        assert graded_profile(windows, unit)
    except AssertionError:
        assert verify_central_profile(windows, unit)
    return degree


def check_start(arguments):
    windows, cuts, start = arguments
    return fitting_support(windows, cuts, start, 0, prime=3, radicand="-1")


def verify():
    windows = integral_windows()
    cuts, bad = path_and_center(windows)
    rules = [[(h, {w: encode(c) for w, c in row.items()}) for h, row in side]
             for side in windows]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for support in pool.map(check_start, [(rules, cuts, p) for p in range(4)]):
            for factor in support:
                bad.update(coefficients(f) for f, _ in R(list(factor)).factor()[1])
        survivors = remaining_orbits(bad)
        degrees = sorted(pool.map(check_point, survivors))
    print("Four length-13 central cuts over F9[s], residual orbit degrees", degrees)


if __name__ == "__main__":
    verify()
