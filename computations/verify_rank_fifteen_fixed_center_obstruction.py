"""Test three central hyperplanes on the quadratic exceptional fiber."""

from functools import cache

from flint import fmpz_mod_poly_ctx, fq_default_ctx
from verify_klein_center_lift import add_multiple, echelon
from verify_rank_fifteen_bivariate_sixth_cut import quotient_profile, sections
from verify_rank_fifteen_bivariate_word_paths import (
    centers0,
    joint_square,
    paths,
    tensor,
)
from verify_rank_fourteen_geometric_reduction import Scalar


def exceptional_point():
    f = fq_default_ctx(modulus=fmpz_mod_poly_ctx(3)([2, 0, 0, 2, 1]), var="a")
    a = f.gen()
    x = a**3 + a * a + 2
    v = a**3 + a * a + a
    c = a * a + 2
    assert x * x + x + 2 == 0 and c * c == x and (v**3 - v) ** 2 == x + 1
    return f, x, c, v


def universal_sections(selected):
    basis, image = paths()
    return sections(basis, centers0(basis, image), selected)


def specialize_sections(f, c, v, symbolic):
    k = 1 / v
    h = c * (k - v)

    def evaluate(poly):
        return sum(
            (int(z) * k**i * h**j for (i, j), z in poly.to_dict().items()), f.zero()
        )

    def specialize(value):
        return Scalar(evaluate(value.n) / evaluate(value.d))

    entries = {w: (e, specialize(z)) for w, (e, z) in tensor.items()}
    joint = specialize(joint_square).value
    assert not joint.is_zero()

    @cache
    def power(exponent):
        assert exponent % 2 == 0
        return Scalar(joint ** (exponent // 2))

    central = [
        [{w: specialize(z) for w, z in row.items()} for row in side]
        for side in symbolic
    ]
    return central, entries, power


def point_certificate():
    f, _, c, v = exceptional_point()
    a = f.gen()
    central, entries, power = specialize_sections(f, c, v, universal_sections(range(5)))
    expected = [4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1, 0]
    ninth = expected[:8] + [23]
    tenth = expected[:9] + [16, 16, 16, 16]
    for label, selected, target in (
        ("fixed", (0, 1, 3, 4), ninth),
        ("exchange", (0, 1, 2, 4), expected),
        ("weighted", (0, 1, 3, 4), tenth),
    ):
        cuts = [[side[i].copy() for i in selected] for side in central]
        if label == "weighted":
            for start, side in enumerate(cuts):
                for i, row in enumerate(side):
                    add_multiple(row, central[start][2], Scalar(a ** (i + 1)))
        for start in range(4):
            profile = quotient_profile(
                cuts, start, len(target), Scalar(f.one()), entries, power, echelon
            )
            rows = [r for n, r in profile]
            if label == "fixed":
                assert rows[-1] == ([5, 6, 6, 6] if start % 2 else [6, 6, 6, 5])
            assert list(map(sum, rows)) == target
            print("Sections", label, "start", start, "dimensions", target, flush=True)


if __name__ == "__main__":
    point_certificate()
