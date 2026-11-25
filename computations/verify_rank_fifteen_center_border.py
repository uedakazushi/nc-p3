"""Quotients in arbitrary central directions from a single five-center reduction."""

from functools import cache

from central_border_paths import central_border
from verify_klein_center_lift import add_multiple, echelon
from verify_rank_fifteen_bivariate_sixth_cut import quotient_profile
from verify_rank_fifteen_fixed_center_obstruction import (
    exceptional_point,
    specialize_sections,
    universal_sections,
)
from verify_rank_fifteen_genus_two_unit_line import tensor_windows
from verify_rank_fifteen_quadratic_parabolic_pencils import beta, word_parity
from verify_rank_fourteen_geometric_reduction import Scalar


def border_data(centers, entries, power, start):
    """Remove word-parity factors before using multiplication over the rational base."""

    def delta(n, h, a):
        return (int((start + n) % 4 == 2 and a == 3) + beta(h, a)) % 2

    def accumulated(word):
        value = h = 0
        for n, a in enumerate(word):
            value += delta(n, h, a)
            h ^= a
        return value

    weighted = [
        {w: z * power(word_parity(start, w) - accumulated(w)) for w, z in form.items()}
        for form in centers[start]
    ]
    windows = tensor_windows(entries)

    @cache
    def relations(n, h):
        rows = []
        for g, relation in windows[(start + n - 2) % 4]:
            k = h ^ g
            row = {}
            parity = None
            for (a, b), (odd, z) in relation.items():
                exponent = odd - delta(n - 2, k, a) - delta(n - 1, k ^ a, b)
                if parity is None:
                    parity = exponent % 2
                assert exponent % 2 == parity
                row[a, b] = z * power(exponent - parity)
            rows.append((g, row))
        return rows

    return relations, weighted


def hyperplane(centers, direction):
    pivot = next(i for i, b in enumerate(direction) if b)
    cuts = []
    for side in centers:
        rows = []
        for j, form in enumerate(side):
            if j != pivot:
                row = form.copy()
                add_multiple(row, side[pivot], -direction[j] / direction[pivot])
                rows.append(row)
        cuts.append(rows)
    return cuts


def degree_bounds():
    full = [
        [1, 0, 0, 0],
        [1] * 4,
        [2, 2, 2, 4],
        [5] * 4,
        [6, 8, 8, 8],
        [9] * 4,
        [10, 10, 10, 4],
        [5] * 4,
    ]
    target = [[5] * 4, [2, 2, 2, 4], [1] * 4, [1, 0, 0, 0], [0] * 4]
    maxima = []
    for n, wanted in enumerate(target, 9):
        row = []
        for h in range(4):
            degrees = [
                (n - i) // 4 for i in range(n % 4, 8, 4) for _ in range(full[i][h])
            ]
            row.append(sum(sorted(degrees, reverse=True)[: len(degrees) - wanted[h]]))
        maxima.append(row)
    assert maxima == [[6] * 4, [12, 12, 12, 8], [14] * 4, [13, 16, 16, 16], [21] * 4]


def verify():
    degree_bounds()
    field, _, c, v = exceptional_point()
    centers, entries, power = specialize_sections(
        field, c, v, universal_sections(range(5))
    )
    one = Scalar(field.one())
    zero = Scalar(field.zero())
    directions = [[one if i == j else zero for i in range(5)] for j in range(5)]
    a = field.gen()
    directions.extend(
        [
            [one, Scalar(a), zero, Scalar(a * a), zero],
            [one, zero, Scalar(a**3), zero, one],
        ]
    )
    for start in range(4):
        relations, weighted = border_data(centers, entries, power, start)
        for direction in directions:
            cuts = hyperplane(centers, direction)
            direct = list(
                quotient_profile(cuts, start, 13, one, entries, power, echelon)
            )
            lifted = list(central_border(relations, weighted, direction, one))
            assert lifted == direct
    print("PASS: homogeneous degree bounds and all 28 complete quotient profiles")


if __name__ == "__main__":
    verify()
