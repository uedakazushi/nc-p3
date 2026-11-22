"""Certify characteristic-zero lifts of every ordinary geometric residue in characteristic 101."""

from concurrent.futures import ProcessPoolExecutor
from functools import cache

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx, nmod_poly
from rank_fourteen_dual import fitting_support, integral_rules
from rank_fourteen_function_field import relations
from verify_central_window_cutoff import quotient_dimensions
from verify_klein_center_lift import central_equations, echelon, field_paths
from verify_rank_fourteen_center import path_chart
from verify_rank_fourteen_central_cut import central_sections
from verify_rank_fourteen_descent import slopes

PROFILE = [1, 4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1, 0]


class Scalar:
    """Supply a zero-sensitive truth value to the existing sparse path arithmetic."""

    __slots__ = ("value",)

    def __init__(self, value):
        self.value = value

    def __bool__(self):
        return not self.value.is_zero()

    def __neg__(self):
        return Scalar(-self.value)

    def __add__(self, other):
        return Scalar(
            self.value + (other.value if isinstance(other, Scalar) else other)
        )

    __radd__ = __add__

    def __mul__(self, other):
        return Scalar(
            self.value * (other.value if isinstance(other, Scalar) else other)
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        return Scalar(
            self.value / (other.value if isinstance(other, Scalar) else other)
        )


@cache
def rational_windows():
    return relations()


def evaluate(poly, t, field):
    value = field.zero()
    for coefficient in reversed(list(poly)):
        value = value * t + int(coefficient)
    return value


def specialize(field, t, i, y):
    def entry(value):
        a, b, c, d = (evaluate(poly, t, field) for poly in value.coordinates)
        return Scalar(
            (a + i * b + y * (c + i * d)) / evaluate(value.denominator, t, field)
        )

    return [
        [(h, {word: entry(c) for word, c in rule.items()}) for h, rule in side]
        for side in rational_windows()
    ]


def test_exception(item):
    poly, imaginary = item
    degree = len(poly) - 1
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(101)(list(poly)), var="a")
    t = field.gen() if degree > 1 else field(-poly[0])
    if (t * t - 1).is_zero() or (t + 3).is_zero():
        return False
    r, lam = slopes(t)
    z = r * r
    if (
        z
        * lam
        * (z * z - 1)
        * (lam * lam * z - 1)
        * (lam * lam - 1)
        * (lam * lam * z * z - 1)
    ).is_zero():
        return False
    cubic = 2 * (t - 1) * (t + 1) * (t + 3)
    if not cubic.is_square():
        field = fq_default_ctx(101, 2 * degree, "a")
        t = fq_default_poly_ctx(field)(list(poly)).roots()[0][0]
        cubic = 2 * (t - 1) * (t + 1) * (t + 3)
    unit = Scalar(field.one())
    windows = specialize(field, t, field(imaginary), cubic.sqrt())
    return verify_central_profile(windows, unit)


def verify_central_profile(windows, unit):
    for side in windows:
        assert (
            len(
                echelon(
                    [
                        {4 * a + b: c for (a, b), c in row.items() if c}
                        for _, row in side
                    ]
                )
            )
            == 6
        )
    basis, image = field_paths(windows, unit)
    assert all(
        sum(len(basis[p, n, h]) for h in range(4)) == expected
        for p in range(4)
        for n, expected in ((4, 35), (5, 56))
    )
    kernels = []
    for h, expected in enumerate((3, 2, 0, 0)):
        rows, size = central_equations(basis, image, h)
        pivots = echelon(rows)
        free = [j for j in range(size) if j not in pivots]
        assert len(free) == expected
        kernels.append((pivots, free))

    def sections(choices):
        out = []
        for h, pairs in enumerate((choices, ((1, 0), (0, 1)))):
            pivots, free = kernels[h]
            for pair in pairs:
                vector = {j: c * unit for j, c in zip(free, pair)}
                vector.update(
                    {
                        j: -sum(
                            (c * row.get(k, 0) for k, c in zip(free, pair)), 0 * unit
                        )
                        for j, row in pivots.items()
                    }
                )
                parts = []
                offset = 0
                for p in range(4):
                    words = basis[p, 4, h]
                    parts.append(
                        {
                            w: vector[offset + k]
                            for k, w in enumerate(words)
                            if vector[offset + k]
                        }
                    )
                    offset += len(words)
                out.append(parts)
        return out

    rules = [[rule for _, rule in side] for side in windows]
    for choices in (
        ((1, 3, 5), (2, 7, 11)),
        ((1, 0, 0), (0, 0, 1)),
        ((1, 2, 3), (4, 1, 5)),
        ((1, 3, 8), (2, 9, 1)),
    ):
        cuts = sections(choices)
        if all(quotient_dimensions(rules, cuts, p, unit) == PROFILE for p in range(4)):
            return True
    raise AssertionError("no common central cut found")


def denominator_factors(rules, prime):
    poles = set()
    for side in rules:
        for _, rule in side:
            for value in rule.values():
                denominator = nmod_poly(list(map(int, value.denominator)), prime)
                assert denominator
                poles.update(tuple(map(int, p)) for p, e in denominator.factor()[1])
    return poles


def run():
    windows = rational_windows()
    basis, image, _ = path_chart(windows)
    cuts = central_sections(basis, image)
    poles = denominator_factors(windows + cuts, 101)
    integral_windows = list(map(integral_rules, windows))
    integral_cuts = list(map(integral_rules, cuts))
    with ProcessPoolExecutor(max_workers=4) as pool:
        jobs = {
            (start, i): pool.submit(
                fitting_support, integral_windows, integral_cuts, start, i
            )
            for i in (10, 91)
            for start in (0, 1)
        }
        exceptions = []
        for i in (10, 91):
            bad = poles.union(*(jobs[start, i].result() for start in (0, 1)))
            exceptions.extend(
                (poly, i) for poly in sorted(bad, key=lambda q: (len(q), q))
            )
            print("Dual Fitting support including poles", i, len(bad), flush=True)
        for item, ordinary in zip(exceptions, pool.map(test_exception, exceptions)):
            print("Residual orbit", *item, "ordinary", ordinary, flush=True)
    print(
        "All ordinary geometric residues in characteristic 101 admit regular characteristic-zero lifts.",
        flush=True,
    )


if __name__ == "__main__":
    run()
