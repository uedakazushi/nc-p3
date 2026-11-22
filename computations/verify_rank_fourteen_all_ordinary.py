"""Saturate the central-cut relations in the two remaining disks over 101."""

from functools import cache
from itertools import product

from flint import ctx, fmpq_poly, fmpq_series
from rank_fourteen_function_field import ONE, relations
from verify_klein_center_lift import add_multiple
from verify_rank_fourteen_center import path_chart
from verify_rank_fourteen_central_cut import central_sections, quotient_layers

INITIAL_PRECISION = 27
DIMENSIONS = (
    (1, 0, 0, 0),
    (1, 1, 1, 1),
    (2, 2, 3, 3),
    (5,) * 4,
    (7, 8, 8, 8),
    (10,) * 4,
    (12, 12, 10, 10),
    (10,) * 4,
    (7, 8, 8, 8),
    (5,) * 4,
    (2, 2, 3, 3),
    (1,) * 4,
    (1, 0, 0, 0),
    (0,) * 4,
)


def gaussian_jets(radicand=-1, initial_precision=INITIAL_PRECISION, residue_prime=101):
    """Truncated quadratic-field power series with a run-local precision."""

    class Jet:
        __slots__ = ("imag", "real")
        precision = initial_precision

        def __init__(self, real, imag=None):
            self.real = real
            self.imag = imag

        @classmethod
        def constant(cls, n):
            return cls(fmpq_poly([n]))

        def __bool__(self):
            return bool(self.real) or bool(self.imag)

        def __neg__(self):
            return Jet(-self.real, -self.imag if self.imag else None)

        def __add__(self, other):
            if isinstance(other, int):
                return Jet(self.real + other, self.imag) if other else self
            imag = (
                self.imag + other.imag
                if self.imag and other.imag
                else self.imag or other.imag
            )
            return Jet(self.real + other.real, imag)

        __radd__ = __add__

        def __sub__(self, other):
            return self + -other

        def __rsub__(self, other):
            return -self + other

        def __mul__(self, other):
            if isinstance(other, int):
                return Jet(self.real * other, self.imag * other if self.imag else None)
            n = self.precision
            real = self.real.mul_low(other.real, n)
            if self.imag and other.imag:
                return Jet(
                    real + radicand * self.imag.mul_low(other.imag, n),
                    self.real.mul_low(other.imag, n) + self.imag.mul_low(other.real, n),
                )
            if self.imag:
                return Jet(real, self.imag.mul_low(other.real, n))
            if other.imag:
                return Jet(real, self.real.mul_low(other.imag, n))
            return Jet(real)

        __rmul__ = __mul__

        def __truediv__(self, other):
            return self * other.inverse()

        def unit(self):
            real = self.real[0]
            imag = self.imag[0] if self.imag else 0
            norm = real * real - radicand * imag * imag
            assert int(norm.denominator) % residue_prime
            return bool(int(norm.numerator) % residue_prime)

        def inverse(self):
            assert self.unit()

            def reciprocal(g):
                return fmpq_poly(list(1 / fmpq_series(g, prec=self.precision)))

            if not self.imag:
                return Jet(reciprocal(self.real))
            n = self.precision
            d = reciprocal(
                self.real.mul_low(self.real, n)
                - radicand * self.imag.mul_low(self.imag, n)
            )
            return Jet(self.real.mul_low(d, n), -self.imag.mul_low(d, n))

        def divide_by_parameter(self):
            assert not self.real[0] and not (self.imag and self.imag[0])
            return Jet(
                self.real.right_shift(1),
                self.imag.right_shift(1) if self.imag else None,
            )

    return Jet


def integral_rows(group):
    common = ONE
    for row in group:
        for value in row.values():
            common *= value.denominator // common.gcd(value.denominator)
    return [
        {
            word: tuple(
                poly * (common // value.denominator) for poly in value.coordinates
            )
            for word, value in row.items()
        }
        for row in group
    ]


def integral_presentation():
    windows = relations()
    basis, image, _ = path_chart(windows)
    cuts = central_sections(basis, image)
    return clear_presentation(windows, cuts)


def clear_presentation(windows, cuts):
    integral_windows = [
        [(h, integral_rows([row])[0]) for h, row in side] for side in windows
    ]
    integral_cuts = [[] for _ in range(4)]
    for j in range(4):
        h = cuts[0][j][0]
        # A single scalar for the four components preserves commutation with every arrow.
        rows = integral_rows([side[j][1] for side in cuts])
        for p, row in enumerate(rows):
            integral_cuts[p].append((h, row))
    return integral_windows, integral_cuts


def expand(presentation, Jet, component, sheet):
    i = Jet(fmpq_poly([]), fmpq_poly([1]))
    u = Jet(fmpq_poly([0, 1]))
    t = -1 + 2 * i + u
    cubic = 2 * (t - 1) * (t + 1) * (t + 3)
    y = sheet * (2 * i) * (2 + 2 * i)
    half = Jet.constant(2).inverse()
    for _ in range(INITIAL_PRECISION.bit_length()):
        y = (y + cubic * y.inverse()) * half
    assert not (y * y - cubic)

    @cache
    def evaluate(coefficients):
        value = Jet.constant(0)
        for coefficient in reversed(coefficients):
            value = value * t + coefficient
        return value

    def entry(coordinates):
        a, b, c, d = (evaluate(tuple(map(int, g))) for g in coordinates)
        return a + component * i * b + y * (c + component * i * d)

    return [
        [
            [
                (h, {word: v for word, coords in row.items() if (v := entry(coords))})
                for h, row in side
            ]
            for side in family
        ]
        for family in presentation
    ]


def saturated_projection(rows, size, target, Jet, *, select_pivot=None):
    pivots = []
    orders = []
    loss = 0
    while len(pivots) < size - target:
        rows = [row for row in rows if row]
        assert rows
        if select_pivot is not None:
            pivot = select_pivot(rows)
        else:
            pivot = next(
                (
                    (j, k)
                    for j in sorted(range(len(rows)), key=lambda j: len(rows[j]))
                    for k in sorted(rows[j])
                    if rows[j][k].unit()
                ),
                None,
            )
        if pivot is None:
            loss += 1
            assert loss < Jet.precision
            rows = [
                {k: v.divide_by_parameter() for k, v in row.items()} for row in rows
            ]
            continue
        orders.append(loss)
        j, k = pivot
        row = rows.pop(j)
        inverse = row.pop(k).inverse()
        normalized = {col: v * inverse for col, v in row.items()}
        pivots.append((k, normalized))
        for other in rows:
            if k in other:
                add_multiple(other, normalized, -other.pop(k))
    used = {j for j, _ in pivots}
    free = [j for j in range(size) if j not in used]
    images = {j: {k: Jet.constant(1)} for k, j in enumerate(free)}
    for j, row in reversed(pivots):
        result = {}
        for k, value in row.items():
            add_multiple(result, images[k], -value)
        images[j] = result
    report = getattr(Jet, "record_pivot_orders", None)
    if report is not None:
        report(tuple(orders))
    return free, [images[j] for j in range(size)], loss


def check_start(presentation, component, sheet, start):
    Jet = gaussian_jets()
    windows, cuts = expand(presentation, Jet, component, sheet)
    losses = []

    def project(n, h, rows, size):
        free, projection, loss = saturated_projection(rows, size, DIMENSIONS[n][h], Jet)
        losses.append(loss)
        return free, projection

    for n in quotient_layers(windows, cuts, start, 13, Jet.constant(1), project):
        # After division each surviving row is correct to this common lower precision.
        Jet.precision -= max(losses)
        assert Jet.precision > 0
        losses.clear()
    return Jet.precision


if __name__ == "__main__":
    ctx.cap = INITIAL_PRECISION
    presentation = integral_presentation()
    for component, sheet, start in product((1, -1), (1, -1), range(4)):
        remaining = check_start(presentation, component, sheet, start)
        print(
            component,
            sheet,
            start,
            "certified through length 13, precision",
            remaining,
            flush=True,
        )
