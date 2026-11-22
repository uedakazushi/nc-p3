"""Degree-six charts and closed-point central cuts on the coordinate ends."""

from itertools import product

from flint import fmpz_poly, nmod

from rank_fourteen_function_field import (
    CurveElement,
    ONE,
    VARIABLE,
    ZERO,
    irreducible_factors,
)
from rank_fifteen_octic import finite_windows
from verify_klein_center_lift import (
    add_multiple,
    central_equations,
    echelon,
    field_paths,
)
from verify_rank_fifteen_quadratic_dual import sixth_insertion
from verify_rank_fourteen_geometric_reduction import verify_central_profile


class QuadraticChart(CurveElement):
    """Use the real and radical coordinates of the existing biquadratic field."""

    precision = 1

    @classmethod
    def constant(cls, value):
        return cls.lift(value)

    def inverse(self):
        a, imaginary, b, mixed = self.coordinates
        assert not imaginary and not mixed
        norm = a * a - self.radicand * b * b
        if norm.degree() > 0:
            self.inverted.add(tuple(map(int, norm)))
        return type(self)(
            (self.denominator * a, ZERO, -self.denominator * b, ZERO), norm
        )


class EllipticChart(QuadraticChart):
    radicand = -(VARIABLE**3) + 14 * VARIABLE**2 - VARIABLE
    inverted = set()


class GaussianChart(QuadraticChart):
    radicand = -ONE
    inverted = set()


def elliptic_frame(t, y):
    return (
        t,
        (4 * t + y) / (t * (t + 1)),
        2 * t**0,
        (t + 1) / (2 * t),
        0 * t,
        4 / (t + 1),
    )


def gaussian_frame(t, i):
    lam = i * t
    z = (1 - t) / (t * (t + 1))
    k = (lam + 1) * (lam * z + 1) / (4 * lam * z)
    m = -(lam - 1) * (lam * z - 1) / (4 * lam * z)
    assert m == i * k
    return z, lam, 1 + i, 2 * k, 1 - i, 0 * t


def sparse_chart(rows, reverse):
    """Choose a shortest row, then its first or last column."""
    pivots = []
    rows = [dict(row) for row in rows if row]
    while rows:
        j = min(range(len(rows)), key=lambda j: len(rows[j]))
        row = rows.pop(j)
        k = (max if reverse else min)(row)
        inverse = row.pop(k).inverse()
        scaled = {c: v * inverse for c, v in row.items()}
        for other in rows:
            if k in other:
                add_multiple(other, scaled, -other.pop(k))
        pivots.append((k, scaled))
        rows = [r for r in rows if r]
    resolved = {}
    for k, row in reversed(pivots):
        for col in list(row):
            if col in resolved:
                add_multiple(row, resolved[col], -row.pop(col))
        resolved[k] = row
    return resolved


def check_paths_and_center(windows, unit, reducer):
    basis, image = field_paths(windows, unit, reducer)
    for p, h in product(range(2), (0, 2)):
        rows, width = sixth_insertion(windows, basis, image, p, h)
        assert (width, len(rows)) == (56, 54 if h == 0 else 51)
        assert len(reducer(rows)) == (36 if h == 0 else 34)

    for h, expected in enumerate((3, 2, 0, 0)):
        rows, size = central_equations(basis, image, h)
        assert size - len(reducer(rows)) == expected


def chart_factors(field, frame, reverse):
    field.inverted.clear()
    x = field.lift(VARIABLE)
    radical = field((ZERO, ZERO, ONE, ZERO))
    unit = field.lift(1)
    windows = finite_windows(frame(x, radical), unit, lambda x: x)
    check_paths_and_center(windows, unit, lambda rows: sparse_chart(rows, reverse))
    return set().union(
        *(irreducible_factors(fmpz_poly(list(poly))) for poly in field.inverted)
    )


# In both sheets of the double cover these parameters are nonordinary or
# outside the affine multiplier chart. Keep ordinary rank-fourteen points.
ELLIPTIC_BOUNDARY = {(0, 1), (1, 1), (-1, 1), (1, -14, 1), (1, -6, 1), (1, 10, 1)}
GAUSSIAN_BOUNDARY = {(0, 1), (1, 1), (-1, 1), (-1, 2, 1), (1, 0, 1), (-1, -2, 1)}
ELLIPTIC_WITNESSES = (
    (103, 92, 7),
    (103, 32, 30),
    (103, 28, 29),
    (103, 29, 5),
    (157, 130, 20),
    (103, 57, 25),
    (101, 96, 50),
    (101, 95, 25),
    (103, 53, 20),
    (101, 65, 14),
    (101, 88, 43),
)
GAUSSIAN_WITNESSES = (
    (109, 60, 33),
    (109, 89, 33),
    (109, 62, 33),
    (109, 58, 33),
    (101, 35, 10),
    (101, 75, 10),
    (101, 99, 10),
)


def verify(field, frame, boundary, sizes, witnesses):
    first, second = (chart_factors(field, frame, reverse) for reverse in (False, True))
    common = first & second
    assert (len(first), len(second), len(common)) == sizes
    residual = sorted(common - boundary, key=lambda g: (len(g), g))
    for g, (p, root, ordinate) in zip(residual, witnesses, strict=True):
        poly = fmpz_poly(list(g))
        assert not int(poly(root)) % p and int(poly.derivative()(root)) % p
        assert ordinate % p and not (ordinate**2 - int(field.radicand(root))) % p
        unit = nmod(1, p)
        for sign in (1, -1):
            values = frame(nmod(root, p), nmod(sign * ordinate, p))
            z, lam = values[:2]
            assert (
                z
                * lam
                * (z * z - 1)
                * (lam * lam * z - 1)
                * (lam * lam - 1)
                * (lam * lam * z * z - 1)
            )
            windows = finite_windows(values, unit, lambda x: x)
            check_paths_and_center(windows, unit, echelon)
            assert verify_central_profile(windows, unit)
    print(
        field.__name__,
        "pivot factors and intersection",
        sizes,
        "certified ordinary factors",
        len(residual),
    )


if __name__ == "__main__":
    verify(
        EllipticChart,
        elliptic_frame,
        ELLIPTIC_BOUNDARY,
        (56, 50, 17),
        ELLIPTIC_WITNESSES,
    )
    verify(
        GaussianChart,
        gaussian_frame,
        GAUSSIAN_BOUNDARY,
        (74, 56, 13),
        GAUSSIAN_WITNESSES,
    )
