"""Functions, relation spaces and central sections on the rank-fifteen ramification curves."""

from functools import cache

from rank_fourteen_function_field import (
    CurveElement,
    ONE,
    VARIABLE,
    ZERO,
    evaluate_relations,
    irreducible_factors,
)
from verify_fixed_geometric_quartic_locus import a, e, f
from verify_klein_center_lift import field_paths
from verify_rank_fifteen_covers import frame_on_branch
from verify_rank_fourteen_all_ordinary import clear_presentation
from verify_rank_fourteen_central_cut import central_sections
from verify_uniform_fixed_quartic_residuals import b
from verify_uniform_fixed_quartic_support import l, z


class BranchElement(CurveElement):
    __slots__ = ()
    radicand = -(VARIABLE**4) + 14 * VARIABLE**2 - 1


T = BranchElement.lift(VARIABLE)
Y = BranchElement((ZERO, ZERO, ONE, ZERO))
UNIT = BranchElement.lift(1)


@cache
def windows(color):
    frame = frame_on_branch(T, Y, color)
    return evaluate_relations(dict(zip((z, l, a, b, e, f), frame)), UNIT)


@cache
def presentation(color):
    rules = windows(color)
    basis, image = field_paths(rules, UNIT)
    assert all(
        sum(len(basis[p, n, h]) for h in range(4)) == d
        for p in range(4)
        for n, d in ((4, 35), (5, 56))
    )
    return clear_presentation(rules, central_sections(basis, image, UNIT))


def boundary_factors():
    z, l, *_ = frame_on_branch(T, Y, 0)
    polys = set()
    for value in (
        T,
        1 + T * T,
        UNIT.lift(BranchElement.radicand),
        z,
        l,
        z * z - 1,
        l * l - 1,
        l * l * z - 1,
        l * l * z * z - 1,
    ):
        polys.update(irreducible_factors(value.norm_numerator()))
        polys.update(irreducible_factors(value.denominator))
    return polys
