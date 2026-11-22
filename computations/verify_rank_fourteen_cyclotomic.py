"""Prove the arithmetic tests for all ordinary cyclotomic values of the rank-14 parameter."""

from concurrent.futures import ProcessPoolExecutor

from flint import fmpz_mod_poly_ctx, fq_default_ctx, nmod_poly
from rank_fourteen_dual import fitting_support, integral_rules
from rank_fourteen_function_field import excluded_factors, relations
from sympy import divisors, n_order, totient
from verify_rank_fourteen_center import path_chart
from verify_rank_fourteen_central_cut import central_sections
from verify_rank_fourteen_geometric_reduction import denominator_factors


def boundary(prime):
    return {
        tuple(map(int, p))
        for f in excluded_factors()
        for p, e in nmod_poly(list(f), prime).factor()[1]
    }


def boundary_orders():
    polys = boundary(101)
    size = sum(len(f) - 1 for f in polys)
    orders = set()
    for f in polys:
        degree = len(f) - 1
        field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(101)(list(f)))
        t = field(-f[0]) if degree == 1 else field.gen()
        orders.add(next(d for d in divisors(101**degree - 1) if (t**d).is_one()))
    assert size == 13 and orders == {1, 2, 25, 51, 100, 3400, 10200}
    assert all(totient(d) > size for d in orders - {1, 2})


def auxiliary_prime():
    assert n_order(109, 101) == 100 and 33**2 % 109 == 108
    windows = relations()
    basis, image, _ = path_chart(windows)
    cuts = central_sections(basis, image)
    bad = boundary(109) | denominator_factors(windows + cuts, 109)
    integral_windows = list(map(integral_rules, windows))
    integral_cuts = list(map(integral_rules, cuts))
    with ProcessPoolExecutor(max_workers=2) as pool:
        jobs = [
            pool.submit(fitting_support, integral_windows, integral_cuts, p, 33, 109)
            for p in (0, 1)
        ]
        for job in jobs:
            bad.update(job.result())
    largest = max(len(f) - 1 for f in bad)
    assert largest < 100
    print("At 109 the largest excluded factor has degree", largest, flush=True)


if __name__ == "__main__":
    boundary_orders()
    auxiliary_prime()
