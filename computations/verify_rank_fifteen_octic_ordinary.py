"""Four-central-form certificates on all ordinary geometric octic residues at 109."""

from concurrent.futures import ProcessPoolExecutor
from functools import cache
from itertools import product
from math import prod

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx
from sympy import sqrt_mod

from rank_fifteen_octic import finite_boundary, finite_windows, octic_frame
from rank_fifteen_octic_modular import F, ONE, P, Function, poly, windows
from verify_klein_center_lift import central_equations, echelon, field_paths
from verify_rank_fourteen_center import rank_minor
from verify_rank_fourteen_central_cut import (
    central_sections,
    exact_projection,
    quotient_layers,
)
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile

INVARIANT_PAIRS = (((1, 0, 0), (0, 0, 1)), ((1, 3, 5), (2, 7, 11)))


@cache
def factors(coefficients):
    assert coefficients
    return frozenset(tuple(map(int, f)) for f, _ in poly(coefficients).factor()[1])


def support(value):
    return factors(tuple(map(int, value)))


@cache
def path_chart(color, j):
    rules = windows(color, j)
    unit = Function(ONE)
    bad = set()
    for side in rules:
        for _, row in side:
            for value in row.values():
                bad.update(support(value.d))

    def reduce(rows):
        rank, det = rank_minor(rows, unit)
        bad.update(support(det.norm_numerator()) | support(det.d))
        pivots = echelon(rows)
        assert rank == len(pivots)
        for row in pivots.values():
            for value in row.values():
                bad.update(support(value.d))
        return pivots

    basis, image = field_paths(rules, unit, reduce)
    assert all(
        sum(len(basis[p, n, h]) for h in range(4)) == d
        for p in range(4)
        for n, d in enumerate((1, 4, 10, 20, 35, 56))
    )
    for h, dimension in enumerate((3, 2, 0, 0)):
        rows, size = central_equations(basis, image, h)
        rank, det = rank_minor(rows, unit)
        assert size - rank == dimension
        bad.update(support(det.norm_numerator()) | support(det.d))
    return basis, image, bad


def cut_polynomial(color, j, start, choice, sample):
    basis, image, bad = path_chart(color, j)
    rules = windows(color, j)
    unit = Function(ONE)
    cuts = central_sections(basis, image, unit, INVARIANT_PAIRS[choice])
    bad = bad.copy()
    for side in cuts:
        for _, row in side:
            for value in row.values():
                bad.update(support(value.d))
    y = int(sqrt_mod(int(F(sample)), P))
    assert y and all(poly(g)(sample) for g in bad)

    def specialize(value):
        return int((value.a(sample) + y * value.b(sample)) / value.d(sample))

    exceptional = prod((poly(g) for g in bad), start=ONE)

    def record(q):
        nonlocal exceptional
        assert q
        for part, _ in q.factor_squarefree()[1]:
            exceptional *= part // part.gcd(exceptional)

    def project(n, character, rows, size):
        det, free, projection = exact_projection(rows, size, unit, specialize, P)
        record(det.norm_numerator())
        record(det.d)
        return free, projection

    for _ in quotient_layers(rules, cuts, start, 13, unit, project):
        pass
    return tuple(map(int, exceptional))


def check_fiber(color, j, factor):
    degree = len(factor) - 1
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(P)(list(factor)), var="a")
    t = field.gen()
    square = t**8 + 14 * t**4 + 1
    assert not square.is_zero()
    if square.is_square():
        y = square.sqrt()
        sheets = (y, -y)
    else:
        field = fq_default_ctx(P, 2 * degree, "a")
        t = fq_default_poly_ctx(field)(list(factor)).roots()[0][0]
        y = (t**8 + 14 * t**4 + 1).sqrt()
        assert y ** (P**degree) == -y
        sheets = (y,)
    for y in sheets:
        frame = octic_frame(t, y, field(j), color)
        assert verify_central_profile(
            finite_windows(frame, field.one(), Scalar), Scalar(field.one())
        )
    return len(sheets)


def run():
    assert F.gcd(F.derivative()) == ONE
    assert all((j**4 - j**2 + 1) % P == 0 for j in (8, 41))
    boundary = set().union(*(support(poly(q)) for q in finite_boundary()))
    with ProcessPoolExecutor(max_workers=3) as pool:
        jobs = {}
        for color, j, start, choice in product(range(2), (8, 41), range(2), range(2)):
            sample = 3 if color == 0 and choice == 1 and start == 0 else 2
            key = color, j, start, choice
            if key == (1, 41, 1, 0):
                sample = 4
            jobs[key] = pool.submit(cut_polynomial, *key, sample)
        for color, j in product(range(2), (8, 41)):
            common = set()
            for start in range(2):
                first, second = (jobs[color, j, start, i].result() for i in range(2))
                common.update(support(poly(first).gcd(poly(second))))
            residual = sorted(
                (q for q in common - boundary if len(q) > 2), key=lambda f: (len(f), f)
            )
            assert (
                len(residual)
                == {(0, 8): 23, (0, 41): 17, (1, 8): 28, (1, 41): 21}[color, j]
            )
            print(
                "Presentation",
                color,
                "j",
                j,
                "higher-degree ordinary factors",
                len(residual),
                flush=True,
            )
            checks = [pool.submit(check_fiber, color, j, factor) for factor in residual]
            for factor, check in zip(residual, checks):
                print(
                    "Verified residue degree",
                    len(factor) - 1,
                    "sheets",
                    check.result(),
                    flush=True,
                )


if __name__ == "__main__":
    run()
