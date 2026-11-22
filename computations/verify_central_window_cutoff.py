"""Degree-16 AS test for rational finite-Frobenius four-windows with marked central quartics."""

from fractions import Fraction
from itertools import combinations, product
from math import comb

from general_cut_tor import CutFrames
from two_window_ext import PathRow, add_scaled, quotient_images
from verify_quartic_obstruction import exact_rank


def apply_word(arrows, degree, vector, word):
    for letter in word:
        result = {}
        for i, c in vector.items():
            add_scaled(result, arrows[degree][4 * i + letter], c)
        degree += 1
        vector = result
    return vector


def central(windows, sections):
    for start in range(4):
        arrows = PathRow(windows, start, 5).arrows
        for section, letter in product(sections, range(4)):
            value = {}
            for word, c in section[start].items():
                add_scaled(value, apply_word(arrows, 0, {0: 1}, word + (letter,)), c)
            for word, c in section[(start + 1) % 4].items():
                add_scaled(value, apply_word(arrows, 0, {0: 1}, (letter,) + word), -c)
            if value:
                return False
    return True


def quotient_dimensions(windows, sections, start):
    dims, arrows = [1], []
    for n in range(1, 14):
        rows = []
        for length in (2, 4):
            if length > n:
                continue
            pos = (start + n - length) % 4
            rules = windows[pos] if length == 2 else [z[pos] for z in sections]
            for i, rule in product(range(dims[n - length]), rules):
                row = {}
                for word, c in rule.items():
                    tail = apply_word(arrows, n - length, {i: 1}, word[:-1])
                    for j, v in tail.items():
                        k = 4 * j + word[-1]
                        row[k] = row.get(k, 0) + c * v
                rows.append({k: v for k, v in row.items() if v})
        size, images = quotient_images(4 * dims[-1], rows)
        dims.append(size)
        arrows.append(images)
    return dims


def koszul_exact(windows, start):
    paths = PathRow(windows, start, 16)
    if paths.dims != [comb(n + 3, 3) for n in range(17)]:
        return False
    frames = CutFrames(windows)
    terms = (1, 4, 6, 4, 1)
    for n in range(1, 17):
        previous_rank = 0
        for j in range(min(n, 4) + 1):
            dim = paths.dims[n - j] * terms[j]
            if j == min(n, 4):
                if previous_rank != dim:
                    return False
                break
            degree = n - j - 1
            split = frames.split("middle", (start + degree) % 4, j + 1, "first")
            columns = []
            for i, by_letter in product(range(paths.dims[degree]), split):
                column = {}
                for letter, parts in enumerate(by_letter):
                    for k, c in paths.arrows[degree][4 * i + letter].items():
                        for b, v in parts.items():
                            idx = k * terms[j] + b
                            column[idx] = column.get(idx, 0) + c * v
                columns.append(column)
            rank = exact_rank(columns)
            if previous_rank + rank != dim:
                return False
            previous_rank = rank
    return True


def certify(windows, sections):
    """Assumes the standard finite Frobenius test. A failed cut leaves AS undecided."""
    windows = [
        [{w: Fraction(c) for w, c in rule.items()} for rule in rs] for rs in windows
    ]
    sections = [
        [{w: Fraction(c) for w, c in rule.items()} for rule in z] for z in sections
    ]
    if len(sections) != 4 or not central(windows, sections):
        return "unmarked"
    expected = [1, 4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1, 0]
    dims = [quotient_dimensions(windows, sections, i) for i in range(4)]
    if any(d[-1] for d in dims):
        return "cut unverified"
    if any(d != expected for d in dims):
        return "non-AS"
    opposite = [
        [{w[::-1]: c for w, c in r.items()} for r in windows[(-i - 2) % 4]]
        for i in range(4)
    ]
    seen = []
    for rules, start in product((windows, opposite), range(4)):
        rotated = rules[start:] + rules[:start]
        if rotated not in seen:
            if not koszul_exact(rotated, 0):
                return "non-AS"
            seen.append(rotated)
    return "AS"


if __name__ == "__main__":
    rules = [{(i, j): 1, (j, i): -1} for i, j in combinations(range(4), 2)]
    windows = [rules] * 4
    sections = [[{(i,) * 4: 1}] * 4 for i in range(4)]
    assert certify(windows, sections) == "AS"
    assert certify(windows, sections[:3] + [sections[0]]) == "cut unverified"
    weights = [[2 + p + 3 * j for j in range(4)] for p in range(4)]
    signed = [
        [
            {
                (i, j): weights[p][i] * weights[(p + 1) % 4][j],
                (j, i): weights[p][j] * weights[(p + 1) % 4][i],
            }
            for i, j in combinations(range(4), 2)
        ]
        for p in range(4)
    ]
    assert certify(signed, sections) == "AS"
    broken = [list(z) for z in sections]
    broken[0][1] = {(0,) * 4: 2}
    assert certify(windows, broken) == "unmarked"
    print(
        "Polynomial four-window certified through 16; failed-cut and failed-centrality controls distinguished."
    )
