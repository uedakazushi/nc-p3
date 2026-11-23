"""Check the four fixed normalized central columns on the exceptional ninth-degree orbit."""

from itertools import product

from verify_klein_center_lift import echelon
from verify_rank_fifteen_bivariate_sixth_cut import quotient_profile
from verify_rank_fifteen_fixed_center_obstruction import (
    exceptional_point,
    specialize_sections,
    universal_sections,
)
from verify_rank_fourteen_geometric_reduction import Scalar


def verify():
    field, x, c0, v0 = exceptional_point()
    symbolic = universal_sections((0, 1, 3, 4))
    for c_sign, v_sign, shift in product((1, -1), (1, -1), (0, 1, -1)):
        c, v = c_sign * c0, v_sign * v0 + shift
        assert c * c == x and (v**3 - v) ** 2 == x + 1
        cuts, entries, power = specialize_sections(field, c, v, symbolic)
        last = []
        for start in range(4):
            profile = quotient_profile(
                cuts, start, 9, Scalar(field.one()), entries, power, echelon
            )
            rows = [row for n, row in profile]
            assert list(map(sum, rows[:8])) == [4, 10, 20, 31, 40, 44, 40, 31]
            assert rows[-1] == ([5, 6, 6, 6] if start % 2 else [6, 6, 6, 5])
            last.append(rows[-1])
        print(
            "Signs",
            c_sign,
            v_sign,
            "shift",
            shift,
            "ninth character ranks",
            last,
            flush=True,
        )


if __name__ == "__main__":
    verify()
