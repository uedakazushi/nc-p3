"""All rank-fifteen branch windows above degree-one ordinary residues at 101."""

from itertools import product

from flint import (
    fmpz_mod_poly_ctx,
    fq_default_ctx,
    fq_default_poly_ctx,
    nmod,
    nmod_poly,
)
from sympy import sqrt_mod

from rank_fifteen_ramification import boundary_factors, windows
from verify_rank_fourteen_geometric_reduction import (
    Scalar,
    evaluate,
    verify_central_profile,
)


def boundary_at_101():
    return {
        tuple(map(int, g))
        for q in boundary_factors()
        for g, _ in nmod_poly(list(q), 101).factor()[1]
    }


def verify_ordinate(color, y, evaluate_poly, wrap, unit):
    def entry(c):
        a, b, d, e = c.coordinates
        assert not b and not e
        return wrap(
            (evaluate_poly(a) + y * evaluate_poly(d)) / evaluate_poly(c.denominator)
        )

    rules = [
        [(h, {w: entry(c) for w, c in row.items()}) for h, row in side]
        for side in windows(color)
    ]
    assert verify_central_profile(rules, unit)


def verify_residue(color, q):
    degree = len(q) - 1
    if degree == 1:
        t = nmod(-q[0], 101)
        ordinates = sqrt_mod(int(-(t**4) + 14 * t * t - 1), 101, all_roots=True)
        if ordinates:
            assert len(ordinates) == 2
            for y in ordinates:
                verify_ordinate(
                    color,
                    nmod(y, 101),
                    lambda p: nmod(int(p(int(t)) % 101), 101),
                    lambda x: x,
                    nmod(1, 101),
                )
            return
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(101)(list(q)), var="a")
    t = field.gen() if degree > 1 else field(-q[0])

    def quartic(t):
        return -(t**4) + 14 * t * t - 1

    value = quartic(t)
    assert not value.is_zero()
    if value.is_square():
        y = value.sqrt()
        ordinates = (y, -y)
    else:
        field = fq_default_ctx(101, 2 * degree, "a")
        t = fq_default_poly_ctx(field)(list(q)).roots()[0][0]
        # Frobenius over the t-residue field interchanges the two ordinates.
        ordinates = (quartic(t).sqrt(),)
    for y in ordinates:
        verify_ordinate(
            color, y, lambda p: evaluate(p, t, field), Scalar, Scalar(field.one())
        )


def representatives():
    boundary = {-q[0] % 101 for q in boundary_at_101() if len(q) == 2}
    assert boundary == {0, 1, 10, 91, 100}
    ordinary = set(range(101)) - boundary
    return sorted(t for t in ordinary if t <= -pow(t, -1, 101) % 101)


def run():
    reps = representatives()
    assert len(reps) == 48
    for color, t in product(range(2), reps):
        verify_residue(color, ((-t) % 101, 1))
        print("Ordinary rational residue", color, t, "certified", flush=True)


if __name__ == "__main__":
    run()
