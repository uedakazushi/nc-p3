"""Universal bivariate path and degree-four center matrices in residue characteristic three."""

from functools import cache, reduce
from itertools import accumulate, product
from math import prod
from operator import xor

import sympy as sp

from flint_rational_function import rational_function_field
from verify_klein_center_lift import add_multiple
from verify_rank_fifteen_genus_two_unit_line import normalize_support
from verify_rank_fifteen_quadratic_parabolic_pencils import word_parity

field = rational_function_field(("K", "H"), prime=3)
K, H = map(field.from_sympy, sp.symbols("K H"))
base = (field.one, (K - 1) / 2, (K + 1) / 2)
squares = (
    (K * K - 1 + 2 * H) / 4,
    (K * K - 1 - 2 * H) / 4,
    (K + 1 + H) / 2,
    (K + 1 - H) / 2,
    (K - 1 + H) / 2,
    (K - 1 - H) / 2,
)
joint_square = prod(squares)
tensor = {}
for word, (sign, powers) in normalize_support().items():
    odd = powers[3] % 2
    coefficient = sign * prod(a**e for a, e in zip(base, powers[:3], strict=True))
    coefficient *= prod(
        a ** ((e - odd) // 2) for a, e in zip(squares, powers[3:], strict=True)
    )
    tensor[word] = odd, coefficient

factors = set()


def record(poly):
    if poly.total_degree() > 0:
        factors.update(str(f) for f, _ in poly.factor()[1])


open_factors = {
    str(f) for value in (K, H, K - 1, K + 1, *squares) for f, _ in value.n.factor()[1]
}
unit_tests = {}


def chart_unit(value):
    key = str(value.n)
    if key not in unit_tests:
        unit_tests[key] = all(str(f) in open_factors for f, _ in value.n.factor()[1])
    return unit_tests[key]


def echelon(rows, *, prefer_chart_units=False):
    rows = [row.copy() for row in rows if row]
    pivot = []
    while rows:

        def score(position):
            row, column = position
            value = rows[row][column]
            return (
                len(rows[row]),
                value.n.total_degree() + value.d.total_degree(),
                column,
            )

        candidates = ((j, k) for j, row in enumerate(rows) for k in row)
        if prefer_chart_units:
            candidates = sorted(candidates, key=score)
            j, k = next(
                ((j, k) for j, k in candidates if chart_unit(rows[j][k])), candidates[0]
            )
        else:
            j, k = min(candidates, key=score)
        row = rows.pop(j)
        entry = row.pop(k)
        record(entry.n)
        inv = field.one / entry
        row = {i: z * inv for i, z in row.items()}
        pivot.append((k, row))
        for r in rows:
            if k in r:
                add_multiple(r, row, -r.pop(k))
        rows = [r for r in rows if r]
    out = {}
    for k, row in reversed(pivot):
        for i in list(row):
            if i in out:
                add_multiple(row, out[i], -row[i])
        row[k] = field.one
        out[k] = row
    return out


@cache
def adjust(ex):
    assert ex % 2 == 0
    return joint_square ** (ex // 2)


def paths():
    basis = {(p, 0, h): [()] if h == 0 else [] for p, h in product(range(4), repeat=2)}
    prepend = {}

    @cache
    def image(p, word):
        if not word:
            return {0: field.one}
        a = word[0]
        tail = word[1:]
        h = reduce(xor, word, 0)
        s = (p + 1) % 4
        projection = prepend[p, len(word), h, a]
        out = {}
        for j, v in image(s, tail).items():
            normal = basis[s, len(tail), h ^ a][j]
            # Convert a reduced tail back to the scaled normal word after prepending.
            e = (
                word_parity(p, word)
                - word_parity(s, tail)
                + word_parity(s, normal)
                - word_parity(p, (a,) + normal)
            )
            add_multiple(out, projection[j], v * adjust(e))
        return out

    for n in range(1, 6):
        for p, h in product(range(4), repeat=2):
            s = (p + 1) % 4
            t = (p + 2) % 4
            words = [(a,) + w for a in range(4) for w in basis[s, n - 1, h ^ a]]
            off = [sum(w[0] < a for w in words) for a in range(5)]
            rows = []
            if n >= 2:
                for g, k in product(range(4), repeat=2):
                    for tail in basis[t, n - 2, h ^ g]:
                        row = {}
                        shift = None
                        for a in range(4):
                            b = a ^ g
                            w = (a, b, k, k ^ g)
                            w = w[-p:] + w[:-p]
                            if w not in tensor:
                                continue
                            exponent, entry = tensor[w]
                            for j, z in image(s, (b,) + tail).items():
                                nr = basis[s, n - 1, h ^ a][j]
                                ex = (
                                    exponent
                                    - word_parity(s, (b,) + tail)
                                    + word_parity(s, nr)
                                    - word_parity(p, (a,) + nr)
                                )
                                if shift is None:
                                    shift = ex % 2
                                col = off[a] + j
                                row[col] = row.get(col, 0) + entry * z * adjust(
                                    ex - shift
                                )
                        rows.append({j: z for j, z in row.items() if z})
            piv = echelon(rows)
            free = [j for j in range(len(words)) if j not in piv]
            ids = {j: i for i, j in enumerate(free)}
            basis[p, n, h] = [words[j] for j in free]
            proj = [
                {ids[j]: field.one}
                if j in ids
                else {ids[k]: -v for k, v in piv[j].items() if k != j}
                for j in range(len(words))
            ]
            for a in range(4):
                prepend[p, n, h, a] = proj[off[a] : off[a + 1]]
        print("Verified all four starts at uncut length", n, flush=True)
        assert all(
            tuple(len(basis[p, n, h]) for h in range(4))
            == (
                (1, 1, 1, 1),
                (2, 2, 2, 4),
                (5, 5, 5, 5),
                (11, 8, 8, 8),
                (14, 14, 14, 14),
            )[n - 1]
            for p in range(4)
        )
    return basis, image


def centers0(basis, image):
    dims = [len(basis[p, 4, 0]) for p in range(4)]
    offsets = [0, *accumulate(dims)]
    rows = []
    for p, a in product(range(4), repeat=2):
        s = (p + 1) % 4
        block = [{} for _ in basis[p, 5, a]]
        for q, right in ((p, True), (s, False)):
            for j, w in enumerate(basis[q, 4, 0]):
                full = w + (a,) if right else (a,) + w
                e = word_parity(q, w) - word_parity(p, full) - word_parity(p, (a,))
                for k, z in image(p, full).items():
                    block[k][offsets[q] + j] = (1 if right else -1) * z * adjust(e)
        rows.extend(block)
    piv = echelon(rows)
    print("Character-zero quartic center rank", offsets[-1] - len(piv), flush=True)
    assert offsets[-1] - len(piv) == 5
    return piv


def other_center(h, basis, image):
    dims = [len(basis[p, 4, h]) for p in range(4)]
    offsets = [0, *accumulate(dims)]
    width = offsets[-1]
    rows = []
    for p, a in product(range(4), repeat=2):
        s = (p + 1) % 4
        block = [[{}, {}] for w in basis[p, 5, a ^ h]]
        for q, sgn in ((p, 1), (s, -1)):
            for j, w in enumerate(basis[q, 4, h]):
                full = w + (a,) if q == p else (a,) + w
                e = word_parity(q, w) - word_parity(p, full)
                odd = e % 2
                for k, z in image(p, full).items():
                    f = sgn * z * adjust(e - odd)
                    col = offsets[q] + j
                    block[k][odd][col] = f
                    block[k][1 - odd][col + width] = joint_square * f if odd else f
        rows.extend(r for pair in block for r in pair)
    piv = echelon(rows, prefer_chart_units=True)
    print("Verified the vanishing center in character", h, flush=True)
    assert len(piv) == 2 * width


def verify():
    basis, image = paths()
    pivots = centers0(basis, image)
    assert len(pivots) == 39 and all(len(basis[p, 4, 0]) == 11 for p in range(4))
    allowed = {
        str(f) for value in (H, K - 1, K + 1, *squares) for f, _ in value.n.factor()[1]
    }
    assert factors == allowed and len(allowed) == 9
    for h in (1, 2, 3):
        other_center(h, basis, image)
    image.cache_clear()
    assert factors == open_factors and len(factors) == 10
    print(
        "All path and center pivot factors belong to the original open chart.",
        flush=True,
    )


if __name__ == "__main__":
    verify()
