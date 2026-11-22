"""Certify every ordinary prime-field multiplier residue at five small primes."""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import product

from flint import fmpz_mod_poly_ctx, fq_default_ctx, nmod, nmod_poly

from rank_fifteen_octic import finite_windows
from verify_klein_center_lift import echelon, field_paths
from verify_rank_fifteen_pencil_linear_residues import graded_profile
from verify_rank_fifteen_quadratic_dual import sixth_insertion
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile

# Number of representatives with smallest root degree 1, 2 or 4.
EXPECTED = {
    7: (0, 0, 1),
    11: (3, 4, 3),
    13: (6, 6, 6),
    17: (14, 15, 10),
    19: (12, 20, 20),
}


def representatives(p):
    points = {
        (z, l)
        for z, l in product(range(1, p), repeat=2)
        if (z * z - 1) * (l * l - 1) * (z * l * l - 1) * (z * z * l * l - 1) % p
    }
    reps = set()
    for z, l in points:
        orbit = {
            (z, l),
            (z, pow(z * l, -1, p)),
            (pow(z, -1, p), pow(l, -1, p)),
            (pow(z, -1, p), z * l % p),
        }
        assert len(orbit) == 4 and orbit <= points
        reps.add(min(orbit))
    assert len(points) == 4 * len(reps)
    return sorted(reps)


def verify_point(args):
    p, z0, l0 = args
    z, l = nmod(z0, p), nmod(l0, p)
    k = (l + 1) * (l * z + 1) / (4 * l * z)
    m = -(l - 1) * (l * z - 1) / (4 * l * z)
    poly = nmod_poly([int(z * m * m), 0, int(1 - z * k * k), 0, 1], p)
    factor = min(
        (g for g, _ in poly.factor()[1]), key=lambda g: (g.degree(), tuple(map(int, g)))
    )
    if factor.degree() == 1:
        n = -factor[0]
        unit = nmod(1, p)
        frame = (z, l, 1 + n, k + m / n, 1 - n, k - m / n)
        windows = finite_windows(frame, unit, lambda x: x)
    else:
        field = fq_default_ctx(
            modulus=fmpz_mod_poly_ctx(p)(list(map(int, factor))), var="a"
        )
        z, l, k, m = (field(int(x)) for x in (z, l, k, m))
        n = field.gen()
        frame = (z, l, 1 + n, k + m / n, 1 - n, k - m / n)
        windows = finite_windows(frame, field.one(), Scalar)
        unit = Scalar(field.one())
    basis, image = field_paths(windows, unit)
    for start in (0, 1):
        for h, target in ((0, 36), (2, 34)):
            rows, width = sixth_insertion(windows, basis, image, start, h)
            assert width == 56 and len(echelon(rows)) == target
    try:
        assert graded_profile(windows, unit)
    except AssertionError:
        assert verify_central_profile(windows, unit)
    return factor.degree()


def verify():
    with ProcessPoolExecutor(max_workers=3) as pool:
        for p, expected in EXPECTED.items():
            counts = Counter(
                pool.map(verify_point, [(p, z, l) for z, l in representatives(p)])
            )
            assert sum(counts.values()) == sum(expected)
            assert tuple(counts[d] for d in (1, 2, 4)) == expected
            print(
                p,
                4 * sum(expected),
                "ordinary multiplier residues certified",
                flush=True,
            )


if __name__ == "__main__":
    verify()
