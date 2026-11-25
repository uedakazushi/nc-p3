"""Degree-thirteen multiplication pencil after a common three-center cut over F81."""

from flint import fq_default_poly_ctx
from verify_klein_center_lift import add_multiple, echelon
from verify_rank_fifteen_bivariate_sixth_cut import quotient_profile
from verify_rank_fifteen_bivariate_word_paths import word_parity
from verify_rank_fifteen_fixed_center_obstruction import (
    exceptional_point,
    specialize_sections,
    universal_sections,
)
from verify_rank_fourteen_geometric_reduction import Scalar


def determinant(rows, field):
    rows = [r[:] for r in rows]
    d = field.one()
    for i in range(len(rows)):
        j = next((j for j in range(i, len(rows)) if not rows[j][i].is_zero()), None)
        if j is None:
            return field.zero()
        if j != i:
            rows[i], rows[j] = rows[j], rows[i]
            d = -d
        pivot = rows[i][i]
        d *= pivot
        for j in range(i + 1, len(rows)):
            ratio = rows[j][i] / pivot
            if not ratio.is_zero():
                for k in range(i + 1, len(rows)):
                    rows[j][k] -= ratio * rows[i][k]
    return d


def pencil(field, left, right):
    ring = fq_default_poly_ctx(field)
    value, product = ring.zero(), ring.one()
    size = len(left)
    for i in range(size + 1):
        t = sum((i // (3**j) % 3 * field.gen() ** j for j in range(4)), field.zero())
        rows = [
            [a + t * b for a, b in zip(x, y, strict=True)]
            for x, y in zip(left, right, strict=True)
        ]
        value += product * ((determinant(rows, field) - value(t)) / product(t))
        product *= ring.gen() - t
    return value


def verify():
    field, _, c, v = exceptional_point()
    centers, entries, power = specialize_sections(
        field, c, v, universal_sections(range(5))
    )
    cuts = [[side[j] for j in (0, 1, 4)] for side in centers]
    x = fq_default_poly_ctx(field).gen()
    b = 2 * field.gen() ** 3 + field.gen() + 2
    expected = x**8 * (x + 1) ** 4 * (x + b) ** 4
    for start in range(4):

        def inspect(n, basis, image):
            if n <= 9:
                # A simultaneous right annihilator would have zero in all four
                # successor character blocks. Remove the common odd parity
                # of each arrow block, retaining the even power in the base field.
                for h in range(4):
                    rows = []
                    for w in basis[n - 1, h]:
                        row = {}
                        for letter in range(4):
                            full = w + (letter,)
                            e = word_parity(start, w) - word_parity(start, full)
                            scale = power(e - e % 2)
                            for j, z in image(full).items():
                                row[letter, j] = scale * z
                        rows.append(row)
                    assert len(echelon(rows)) == len(rows)
            if n != 13:
                return
            for h in range(4):
                words = basis[9, h]
                width = len(basis[13, h])
                assert width == len(words)
                matrices = []
                for z in (centers[start][3], centers[start][2]):
                    rows = []
                    for tail in words:
                        row = {}
                        for w, coefficient in z.items():
                            e = (
                                word_parity(start, w)
                                + word_parity(start, tail)
                                - word_parity(start, w + tail)
                            )
                            add_multiple(row, image(w + tail), coefficient * power(e))
                        rows.append(
                            [
                                row[j].value if j in row else field.zero()
                                for j in range(width)
                            ]
                        )
                    matrices.append(rows)
                p = pencil(field, *matrices)
                assert width == 16 and p.monic() == expected
                print(
                    "Start", start, "character", h, "determinant verified", flush=True
                )

        profile = list(
            quotient_profile(
                cuts,
                start,
                13,
                Scalar(field.one()),
                entries,
                power,
                echelon,
                inspect_degree=inspect,
            )
        )
        assert [sum(ds) for _, ds in profile] == [
            4,
            10,
            20,
            32,
            44,
            54,
            60,
            63,
            64,
            64,
            64,
            64,
            64,
        ]

    # Independent calculation in the actual quotients at each exceptional root,
    # at a regular finite parameter, and after all five cuts.
    initial = [4, 10, 20, 31, 40, 44, 40, 31]
    for slope, target in (
        (field.zero(), initial + [23, 22, 22, 22, 22]),
        (-field.one(), initial + [20, 16, 16, 16, 16]),
        (-b, initial + [20, 16, 16, 16, 16]),
        (field.one(), initial + [20, 10, 4, 1, 0]),
        (None, [4, 10, 20, 30, 36, 34, 20, 0]),
    ):
        sliced = centers
        if slope is not None:
            sliced = [
                triple + [side[3].copy()]
                for triple, side in zip(cuts, centers, strict=True)
            ]
            for rows, side in zip(sliced, centers, strict=True):
                add_multiple(rows[-1], side[2], Scalar(slope))
        for start in range(4):
            profile = quotient_profile(
                sliced, start, len(target), Scalar(field.one()), entries, power, echelon
            )
            assert [sum(row) for _, row in profile] == target
        print("Direct quotient", slope, target, flush=True)


if __name__ == "__main__":
    verify()
