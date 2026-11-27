"""Certify the four-central-form quotient through length six on the ordinary rank-14 locus."""

from functools import cache, reduce
from operator import xor

from rank_fourteen_function_field import (
    CurveElement,
    excluded_factors,
    irreducible_factors,
    one,
    relations,
)
from verify_klein_center_lift import add_multiple, central_equations, echelon
from verify_rank_fourteen_center import certify_residue, path_chart


def central_sections(basis, image, unit=one, invariant_pairs=((1, 0, 0), (0, 0, 1))):
    sections = [[] for _ in range(4)]
    for h, selected in enumerate((invariant_pairs, ((1, 0), (0, 1)))):
        rows, size = central_equations(basis, image, h)
        pivots = echelon(rows)
        free = [k for k in range(size) if k not in pivots]
        assert len(free) == 3 - h
        for weights in selected:
            assert len(free) == len(weights)
            vector = {k: value * unit for k, value in zip(free, weights) if value}
            for k, row in pivots.items():
                value = -sum(
                    (w * row[col] for col, w in zip(free, weights) if w and col in row),
                    0 * unit,
                )
                if value:
                    vector[k] = value
            offset = 0
            for p, components in enumerate(sections):
                words = basis[p, 4, h]
                components.append(
                    (
                        h,
                        {
                            word: vector[offset + j]
                            for j, word in enumerate(words)
                            if offset + j in vector
                        },
                    )
                )
                offset += len(words)
            assert offset == size
    return sections


def residue(value):
    value = CurveElement.lift(value)
    a, b, c, d = (int(poly(2)) % 101 for poly in value.coordinates)
    return (
        (a + 10 * b + 38 * (c + 10 * d)) * pow(int(value.denominator(2)), -1, 101) % 101
    )


def independent_rows(rows, specialize=residue, prime=101):
    """Choose independent rows at a split prime before computing the determinant in the original field."""
    pivots = {}
    for index in sorted(range(len(rows)), key=lambda j: (len(rows[j]), j)):
        row = {j: v for j, value in rows[index].items() if (v := specialize(value))}
        for j in sorted(pivots):
            if j in row:
                factor = row[j]
                for k, value in pivots[j].items():
                    result = (row.get(k, 0) - factor * value) % prime
                    if result:
                        row[k] = result
                    else:
                        row.pop(k, None)
        if row:
            j = min(row)
            inverse = pow(row[j], -1, prime)
            pivots[j] = {k: value * inverse % prime for k, value in row.items()}
            yield rows[index].copy()


def exact_projection(rows, size, unit=one, specialize=residue, prime=101):
    pivots = {}
    determinant = unit
    for row in independent_rows(rows, specialize, prime):
        for j in sorted(pivots):
            if j in row:
                add_multiple(row, pivots[j], -row[j])
        assert row
        j = min(row)
        determinant *= row[j]
        inverse = row[j].inverse()
        pivots[j] = {k: value * inverse for k, value in row.items()}
    for j in sorted(pivots, reverse=True):
        for k in list(pivots[j]):
            if k > j and k in pivots:
                add_multiple(pivots[j], pivots[k], -pivots[j][k])
    free = [j for j in range(size) if j not in pivots]
    position = {j: k for k, j in enumerate(free)}
    projection = [
        {position[j]: unit}
        if j in position
        else {position[k]: -v for k, v in pivots[j].items() if k != j}
        for j in range(size)
    ]
    return determinant, free, projection


FOUR_QUARTIC_PROFILE = (1, 4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1, 0)


def quotient_layers(windows, sections, start, limit, unit, project, profile=FOUR_QUARTIC_PROFILE):
    """Append one arrow and impose terminal quadratic relations and the central cut."""
    paths = {(0, h): [()] if h == 0 else [] for h in range(4)}
    arrows = {}

    @cache
    def image(word):
        if not word:
            return {0: unit}
        prefix = word[:-1]
        h = reduce(xor, prefix, 0)
        row = {}
        for j, value in image(prefix).items():
            add_multiple(row, arrows[len(prefix), h, word[-1]][j], value)
        return row

    for n in range(1, limit + 1):
        for h in range(4):
            words = [prefix + (a,) for a in range(4) for prefix in paths[n - 1, h ^ a]]
            offsets = [sum(w[-1] < a for w in words) for a in range(5)]
            rows = []
            if n >= 2:
                for g, rule in windows[(start + n - 2) % 4]:
                    for j in range(len(paths[n - 2, h ^ g])):
                        row = {}
                        for (a, b), c in rule.items():
                            for k, v in arrows[n - 2, h ^ g, a][j].items():
                                col = offsets[b] + k
                                row[col] = row.get(col, 0) + c * v
                        rows.append({k: v for k, v in row.items() if v})
            if n == 4:
                for g, section in sections[start]:
                    if g != h:
                        continue
                    row = {}
                    for word, c in section.items():
                        for k, v in image(word[:-1]).items():
                            col = offsets[word[-1]] + k
                            row[col] = row.get(col, 0) + c * v
                    rows.append({k: v for k, v in row.items() if v})
            free, projection = project(n, h, rows, len(words))
            paths[n, h] = [words[j] for j in free]
            for a in range(4):
                arrows[n - 1, h ^ a, a] = projection[offsets[a] : offsets[a + 1]]
        assert sum(len(paths[n, h]) for h in range(4)) == profile[n], (start, n)
        yield n


def quotient_from_start(windows, sections, start, limit=6):
    factors = set()

    def project(n, h, rows, size):
        det, free, projection = exact_projection(rows, size)
        factors.update(irreducible_factors(det.norm_numerator()))
        factors.update(irreducible_factors(det.denominator))
        return free, projection

    for n in quotient_layers(windows, sections, start, limit, one, project):
        print("Verified start", start, "length", n, flush=True)
    return factors


WITNESSES = (
    (8, 157, 141),
    (8, 101, 89),
    (8, 113, 101),
    (14, 101, 51),
    (16, 113, 66),
    (22, 101, 16),
    (16, 109, 94),
    (20, 101, 45),
    (24, 109, 61),
    (24, 109, 90),
    (34, 101, 32),
    (34, 101, 85),
    (52, 109, 90),
    (52, 101, 94),
    (76, 101, 0),
    (106, 101, 17),
    (100, 101, 94),
    (112, 101, 3),
    (108, 101, 86),
    (116, 101, 0),
    (114, 101, 84),
    (128, 101, 56),
)


if __name__ == "__main__":
    windows = relations()
    basis, image, factors = path_chart(windows)
    sections = central_sections(basis, image)
    for rules in windows + sections:
        for _, rule in rules:
            for value in rule.values():
                factors.update(irreducible_factors(value.denominator))
    for start in (0, 1):
        factors.update(quotient_from_start(windows, sections, start))
    factors -= excluded_factors()
    assert len(factors) == len(WITNESSES)
    for polynomial, (degree, prime, root) in zip(sorted(factors), WITNESSES):
        assert len(polynomial) - 1 == degree
        certify_residue(list(polynomial), prime, root)
    print(
        "Uniform central cut through length six; 4768 geometric windows certified in all degrees."
    )
