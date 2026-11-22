"""Uniform center and length-six central cut on the two ordinary ramification curves."""

from functools import cache

from flint import fmpz_poly, nmod, nmod_poly
from sympy import sqrt_mod

from rank_fifteen_ramification import UNIT, boundary_factors, windows
from rank_fourteen_function_field import irreducible_factors
from verify_klein_center_lift import central_equations, echelon, field_paths
from verify_rank_fifteen_covers import frame_on_branch
from verify_rank_fourteen_center import rank_minor
from verify_rank_fourteen_central_cut import (
    central_sections,
    exact_projection,
    quotient_layers,
)
from verify_rank_fourteen_geometric_reduction import verify_central_profile


def shifted_factor(polynomial):
    # Homogenize t -> -1/t, then restore the positive leading coefficient.
    result = tuple((-1) ** j * c for j, c in enumerate(reversed(polynomial)))
    return result if result[-1] > 0 else tuple(-c for c in result)


def exceptional_support(color):
    rules = windows(color)
    factors = set()

    @cache
    def record(coefficients):
        factors.update(irreducible_factors(fmpz_poly(list(coefficients))))

    def store(poly):
        record(tuple(map(int, poly)))

    def minor(det):
        store(det.norm_numerator())
        store(det.denominator)

    def paths(rows):
        rank, det = rank_minor(rows, UNIT)
        minor(det)
        pivots = echelon(rows)
        assert len(pivots) == rank
        for row in pivots.values():
            for value in row.values():
                store(value.denominator)
        return pivots

    basis, image = field_paths(rules, UNIT, paths)
    assert all(
        sum(len(basis[p, n, h]) for h in range(4)) == dimension
        for p in range(4)
        for n, dimension in enumerate((1, 4, 10, 20, 35, 56))
    )
    for h, dimension in enumerate((3, 2, 0, 0)):
        rows, width = central_equations(basis, image, h)
        rank, det = rank_minor(rows, UNIT)
        assert width - rank == dimension
        minor(det)
    cuts = central_sections(basis, image, UNIT)
    for side in rules + cuts:
        for _, row in side:
            for value in row.values():
                store(value.denominator)

    def residue(value):
        a, b, c, d = value.coordinates
        assert not b and not d
        # (t,y)=(4,13) on the quartic at 101 selects rows only.
        return (
            int((a(4) + 13 * c(4)) % 101)
            * pow(int(value.denominator(4)), -1, 101)
            % 101
        )

    def project(n, h, rows, size):
        det, free, projection = exact_projection(rows, size, UNIT, residue)
        minor(det)
        return free, projection

    for start in (0, 1):
        for _ in quotient_layers(rules, cuts, start, 6, UNIT, project):
            pass
    factors -= boundary_factors()
    return sorted(factors | {shifted_factor(f) for f in factors})


@cache
def verify_both_sheets(color, prime, root):
    t = nmod(root, prime)
    roots = sqrt_mod(int(-(t**4) + 14 * t * t - 1), prime, all_roots=True)
    assert len(roots) == 2
    for y0 in roots:
        y = nmod(y0, prime)
        z, lam, *_ = frame_on_branch(t, y, color)
        assert (
            z
            * lam
            * (z * z - 1)
            * (lam * lam - 1)
            * (lam * lam * z - 1)
            * (lam * lam * z * z - 1)
        )

        def entry(value):
            a, b, c, d = value.coordinates
            assert not b and not d
            return (int(a(root) % prime) + y * int(c(root) % prime)) / int(
                value.denominator(root) % prime
            )

        rules = [
            [(h, {word: entry(v) for word, v in row.items()}) for h, row in side]
            for side in windows(color)
        ]
        assert verify_central_profile(rules, nmod(1, prime))


WITNESSES = (
    (
        (101, 14),
        (109, 96),
        (109, 92),
        (101, 89),
        (101, 36),
        (109, 42),
        (101, 93),
        (101, 72),
        (101, 54),
        (109, 77),
        (103, 91),
        (107, 99),
        (103, 43),
        (101, 13),
        (107, 41),
        (107, 88),
        (109, 63),
        (157, 106),
        (101, 87),
        (101, 63),
        (103, 5),
        (103, 73),
        (101, 86),
        (101, 90),
        (101, 87),
        (101, 76),
        (103, 91),
        (109, 64),
        (157, 117),
        (101, 88),
        (127, 120),
        (101, 48),
        (101, 31),
        (107, 62),
        (107, 67),
        (103, 41),
        (107, 71),
        (103, 43),
        (103, 79),
        (103, 19),
        (101, 93),
        (101, 65),
        (107, 60),
        (103, 14),
        (101, 38),
        (101, 97),
        (103, 65),
        (107, 41),
        (103, 22),
    ),
    (
        (101, 60),
        (101, 70),
        (101, 89),
        (101, 69),
        (101, 88),
        (101, 74),
        (101, 72),
        (113, 107),
        (127, 89),
        (101, 87),
        (101, 76),
        (101, 97),
        (107, 88),
        (101, 94),
        (157, 106),
        (101, 63),
        (101, 87),
        (107, 53),
        (101, 86),
        (101, 90),
        (103, 22),
        (107, 60),
        (101, 59),
        (109, 96),
        (101, 29),
        (157, 117),
        (101, 42),
        (101, 88),
        (127, 123),
        (107, 104),
        (101, 88),
        (107, 62),
        (101, 65),
        (113, 19),
        (107, 2),
        (101, 64),
        (103, 98),
        (103, 88),
        (101, 72),
        (101, 93),
        (103, 62),
        (101, 89),
        (109, 42),
        (101, 71),
        (101, 89),
        (103, 55),
        (101, 90),
    ),
)


def run():
    for color, count, points in ((0, 49, 2720), (1, 47, 2616)):
        support = exceptional_support(color)
        assert len(support) == count and 2 * sum(len(f) - 1 for f in support) == points
        assert len(WITNESSES[color]) == count
        for poly, (prime, root) in zip(support, WITNESSES[color]):
            g = nmod_poly(list(poly), prime)
            assert g.degree() == len(poly) - 1 and not g(root) and g.derivative()(root)
            verify_both_sheets(color, prime, root)
        print(
            "Ramification curve",
            color,
            "center and cut through six everywhere;",
            points,
            "regular geometric points.",
            flush=True,
        )


if __name__ == "__main__":
    run()
