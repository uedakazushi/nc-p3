"""Four-center quotient certificates for the complete coordinate-end curves."""

from concurrent.futures import ProcessPoolExecutor
from functools import partial

from flint import fmpz_poly, nmod

from rank_fifteen_octic import finite_windows
from verify_rank_fifteen_coordinate_ends import (
    ELLIPTIC_BOUNDARY,
    EllipticChart,
    ONE,
    VARIABLE,
    ZERO,
    elliptic_frame,
    irreducible_factors,
)
from verify_rank_fifteen_pencil_trace_disks import local_central_sections
from verify_rank_fourteen_all_ordinary import DIMENSIONS, saturated_projection
from verify_rank_fourteen_central_cut import quotient_layers
from verify_rank_fourteen_geometric_reduction import verify_central_profile


def lowest_degree_pivot(rows):
    """Prefer sparse rows and low-degree coefficients in the quadratic field."""

    def score(item):
        j, k = item
        value = rows[j][k]
        degree = max(value.coordinates[i].degree() for i in (0, 2))
        denominator = value.denominator.degree()
        return (
            len(rows[j]) // 4,
            max(degree, denominator),
            len(rows[j]),
            min(degree, denominator),
            k,
        )

    return min(((j, k) for j, row in enumerate(rows) for k in row), key=score)


def central_cut_factors(field, frame):
    field.precision = 1
    field.inverted.clear()
    unit = field.constant(1)
    relations = finite_windows(
        frame(field.lift(VARIABLE), field((ZERO, ZERO, ONE, ZERO))),
        unit,
        lambda value: value,
    )
    project = partial(saturated_projection, select_pivot=lowest_degree_pivot)
    cuts = local_central_sections(
        relations, field, project=project, invariant_pairs=((1, 0, 0), (0, 0, 1))
    )

    def project_paths(degree, h, rows, size):
        free, images, loss = project(rows, size, DIMENSIONS[degree][h], field)
        assert loss == 0
        return free, images

    for start in (0, 1):
        for _ in quotient_layers(relations, cuts, start, 13, unit, project_paths):
            pass
        print(field.__name__, "start", start, "certified through degree 13", flush=True)
    return set().union(
        *(irreducible_factors(fmpz_poly(list(poly))) for poly in field.inverted)
    )


def reciprocal(poly, sign):
    coeffs = list(reversed([v * sign**j for j, v in enumerate(poly)]))
    return tuple(v if coeffs[-1] > 0 else -v for v in coeffs)


def residual_orbits(factors, sign, boundary):
    common = {poly for poly in factors if poly[0] and reciprocal(poly, sign) in factors}
    residual = sorted(
        (g for g in common - boundary if g <= reciprocal(g, sign)),
        key=lambda g: (len(g), g),
    )
    return common, residual


def check_residue(field, frame, poly, witness):
    p, t, y = witness
    f = fmpz_poly(list(poly))
    assert not int(f(t)) % p and int(f.derivative()(t)) % p
    assert y % p and not (y * y - int(field.radicand(t))) % p
    unit = nmod(1, p)
    for sign in (1, -1):
        values = frame(nmod(t, p), nmod(sign * y, p))
        z, lam = values[:2]
        assert (
            z
            * lam
            * (z * z - 1)
            * (lam * lam * z - 1)
            * (lam * lam - 1)
            * (lam * lam * z * z - 1)
        )
        assert verify_central_profile(finite_windows(values, unit, lambda v: v), unit)


ELLIPTIC_WITNESSES = (
    (109, 84, 23),
    (101, 83, 36),
    (107, 96, 19),
    (107, 51, 36),
    (103, 92, 7),
    (157, 130, 20),
    (101, 64, 41),
    (103, 57, 25),
    (101, 96, 50),
    (107, 104, 7),
    (101, 95, 25),
    (101, 69, 26),
    (103, 53, 20),
    (109, 74, 1),
    (103, 78, 35),
    (127, 122, 37),
    (131, 102, 6),
    (181, 140, 27),
    (103, 53, 20),
    (101, 89, 25),
    (109, 80, 2),
    (103, 86, 18),
    (101, 86, 50),
    (131, 79, 15),
    (131, 108, 48),
    (107, 72, 51),
    (103, 6, 30),
    (107, 91, 45),
    (103, 94, 11),
    (127, 105, 38),
    (103, 80, 51),
    (107, 47, 43),
    (113, 70, 32),
    (101, 34, 50),
    (107, 97, 22),
    (113, 52, 34),
    (107, 104, 7),
    (109, 41, 5),
    (103, 38, 29),
    (107, 58, 11),
    (107, 99, 5),
    (109, 50, 23),
    (107, 72, 51),
    (103, 62, 38),
    (103, 92, 7),
    (109, 81, 16),
    (101, 89, 25),
    (107, 96, 19),
    (103, 92, 7),
    (101, 92, 16),
    (101, 92, 16),
    (101, 92, 16),
    (103, 78, 35),
    (107, 51, 36),
    (127, 100, 43),
)


def certify(family, frame, boundary, sign, expected, witnesses):
    factors = central_cut_factors(family, frame)
    confirm_residues(factors, family, frame, boundary, sign, expected, witnesses)


def confirm_residues(factors, family, frame, boundary, sign, expected, witnesses):
    common, residual = residual_orbits(factors, sign, boundary)
    assert (len(factors), len(common), len(residual)) == expected
    check = partial(check_residue, family, frame)
    with ProcessPoolExecutor(max_workers=3) as pool:
        jobs = [
            pool.submit(check, g, w) for g, w in zip(residual, witnesses, strict=True)
        ]
        for job in jobs:
            job.result()
    print(
        family.__name__,
        "pivot factors, common factors, ordinary orbits",
        expected,
        flush=True,
    )


if __name__ == "__main__":
    certify(
        EllipticChart,
        elliptic_frame,
        ELLIPTIC_BOUNDARY,
        1,
        (464, 69, 55),
        ELLIPTIC_WITNESSES,
    )
