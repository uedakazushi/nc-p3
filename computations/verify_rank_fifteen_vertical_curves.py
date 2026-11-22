"""At least one regular infinity end on every smooth rank-fifteen level."""

from concurrent.futures import ProcessPoolExecutor
from functools import cache

import sympy as sp
from flint import fmpz_poly, nmod

from verify_rank_fifteen_pencil_infinity import v, weighted_window
from verify_rank_fifteen_pencil_trace_disks import local_central_sections
from verify_rank_fifteen_relative_trace_centers import Coefficient
from verify_rank_fourteen_all_ordinary import DIMENSIONS, saturated_projection
from verify_rank_fourteen_central_cut import quotient_layers
from verify_rank_fourteen_geometric_reduction import verify_central_profile


def normalized_end(color):
    """Compute periodic arrow scalings from the exponent matrix, allowing roots."""
    windows, units = weighted_window(color, 0, 0)
    generators = (v, v - 1, v + 1, v * v + 1)

    @cache
    def factor(value):
        powers = {g: 0 for g in generators}
        scalar = sp.Integer(1)
        for polynomial, exponent in zip(sp.fraction(sp.cancel(value)), (1, -1)):
            leading, factors = sp.factor_list(polynomial, v)
            scalar *= leading**exponent
            for g, multiplicity in factors:
                powers[g] += exponent * multiplicity
        a, b = map(int, sp.fraction(abs(scalar)))
        assert a > 0 and b > 0 and a & (a - 1) == 0 and b & (b - 1) == 0
        return int(sp.sign(scalar)), [*powers.values(), a.bit_length() - b.bit_length()]

    for value in units:
        factor(value)

    incidence, exponents, signs, inputs = [], [], [], []
    for p, side in enumerate(windows):
        for row_id, (_, row) in enumerate(side):
            (a, b), first = next(iter(row.items()))
            for (i, j), value in row.items():
                line = [0] * 16
                line[4 * p + a] += 1
                line[4 * p + i] -= 1
                nxt = 4 * ((p + 1) % 4)
                line[nxt + b] += 1
                line[nxt + j] -= 1
                incidence.append(line)
                sign, powers = factor(value / first)
                exponents.append(powers)
                signs.append(sign)
                inputs.append((p, row_id, (i, j)))
    matrix, exponents = sp.Matrix(incidence), sp.Matrix(exponents)
    indices = list(matrix.T.rref()[1])
    solution, parameters = matrix[indices, :].gauss_jordan_solve(exponents[indices, :])
    solution = solution.xreplace(dict.fromkeys(parameters, 0))
    rest = exponents - matrix * solution
    assert rest[:, 1:] == sp.zeros(rest.rows, 4)
    out = [[(h, {}) for h, _ in side] for side in windows]
    for (p, row_id, word), sign, power in zip(inputs, signs, rest[:, 0]):
        assert power.is_Integer and power % 2 == 0
        out[p][row_id][1][word] = (sign, int(power // 2))
    return out


class TrackedFraction:
    """Elements of Q(y) recording every inverted numerator on the chosen chart."""

    precision = 1
    inverted = set()

    def __init__(self, value):
        self.value = value

    @classmethod
    def constant(cls, value):
        return cls(Coefficient.constant(value))

    def __bool__(self):
        return bool(self.value)

    def __neg__(self):
        return TrackedFraction(-self.value)

    def __add__(self, other):
        if not isinstance(other, TrackedFraction):
            other = self.constant(other)
        return TrackedFraction(self.value + other.value)

    __radd__ = __add__

    def __mul__(self, other):
        if not isinstance(other, TrackedFraction):
            other = self.constant(other)
        return TrackedFraction(self.value * other.value)

    __rmul__ = __mul__

    def unit(self):
        return bool(self)

    def inverse(self):
        numerator = self.value.num
        if numerator.degree() > 0:
            self.inverted.add(tuple(map(int, numerator)))
        return TrackedFraction(self.value.inverse())

    def __truediv__(self, other):
        return self * other.inverse()


def evaluate(window, monomial):
    return [
        [(h, {w: monomial(*entry) for w, entry in row.items()}) for h, row in side]
        for side in window
    ]


def pivot_factors(window):
    field = TrackedFraction
    field.inverted.clear()
    y = fmpz_poly([0, 1])

    @cache
    def monomial(sign, power):
        return field(Coefficient(sign * y ** max(0, power), y ** max(0, -power)))

    relations = evaluate(window, monomial)
    cuts = local_central_sections(relations, field)
    for start in range(4):

        def project(n, h, rows, size):
            free, images, loss = saturated_projection(
                rows, size, DIMENSIONS[n][h], field
            )
            assert loss == 0
            return free, images

        for _ in quotient_layers(
            relations, cuts, start, 13, field.constant(1), project
        ):
            pass
    return {(0, 1)} | {
        tuple(map(int, factor))
        for poly in field.inverted
        for factor, _ in fmpz_poly(list(poly)).factor()[1]
    }


def reciprocal(poly):
    sign = 1 if poly[0] > 0 else -1
    return tuple(sign * a for a in reversed(poly))


WITNESSES = (
    (
        (101, 3),
        (101, 2),
        (101, 62),
        (107, 91),
        (101, 80),
        (101, 91),
        (103, 56),
        (101, 50),
        (107, 100),
        (107, 46),
        (101, 75),
    ),
    (
        (101, 3),
        (101, 2),
        (101, 71),
        (101, 85),
        (103, 66),
        (101, 6),
        (101, 89),
        (101, 62),
        (101, 36),
        (101, 88),
        (107, 91),
        (103, 86),
        (101, 80),
        (101, 46),
        (107, 66),
        (103, 57),
        (101, 91),
        (101, 67),
        (107, 87),
    ),
)


def certify(color):
    window = normalized_end(color)
    factors = pivot_factors(window)
    common = {f for f in factors if f[0] and reciprocal(f) in factors}
    assert (len(factors), len(common)) == ((994, 17), (923, 33))[color]
    y = fmpz_poly([0, 1])
    # y=0,±1 and xi=1 or xi²+xi+1=0 are outside the smooth genus-five levels.
    excluded = (
        (y * y - 1)
        * (y * y - 6 * y + 1)
        * ((y - 1) ** 4 + 4 * y * (y - 1) ** 2 + 16 * y * y)
    )
    residual = sorted(
        f for f in common if f <= reciprocal(f) and excluded % fmpz_poly(list(f))
    )
    for f, (prime, root) in zip(residual, WITNESSES[color], strict=True):
        poly = fmpz_poly(list(f))
        assert not int(poly(root)) % prime
        assert int(poly.derivative()(root)) % prime and int(excluded(root)) % prime
        value = nmod(root, prime)
        relations = evaluate(window, lambda sign, power: sign * value**power)
        assert verify_central_profile(relations, nmod(1, prime))
    print(
        "Presentation",
        color,
        "pivot factors",
        len(factors),
        "reciprocal support",
        len(common),
        "smooth residual orbits",
        len(residual),
        flush=True,
    )


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=2) as pool:
        list(pool.map(certify, range(2)))
