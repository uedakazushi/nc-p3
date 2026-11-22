"""Certify ordinary characteristic-zero lifts of every nonordinary disk over 101."""

from functools import cache
from itertools import product

from flint import (
    ctx,
    fmpq,
    fmpq_poly,
    fmpz_mod_poly_ctx,
    fq_default_ctx,
    fq_default_poly_ctx,
    nmod_poly,
)

from rank_fifteen_ramification import presentation
from verify_klein_center_lift import echelon
from verify_rank_fourteen_all_ordinary import (
    DIMENSIONS,
    gaussian_jets,
    saturated_projection,
)
from verify_rank_fourteen_central_cut import quotient_layers
from verify_rank_fourteen_geometric_reduction import (
    Scalar,
    evaluate,
    verify_central_profile,
)

# Coordinates are the constant terms (a+b sqrt(d), c+e sqrt(d)).
# At infinity the coordinates are x=1/t and y'=x² y.
CENTERS = {
    "1": (3, 1, 0, 0, 2),
    "0": (-1, 0, 0, 0, 1),
    "-1": (3, -1, 0, 0, 2),
    "pole+": (2, 1, 1, 4, 2),
    "pole-": (2, -1, 1, 4, -2),
    "ram+": (3, 2, 1, 0, 0),
    "i": (-1, 0, 1, 0, 4),
    "ram-": (3, -2, 1, 0, 0),
    "infinity": (-1, 0, 0, 0, 1),
}
PRECISION = {"1": 6, "0": 16, "-1": 6, "ram+": 14, "i": 17, "ram-": 15, "infinity": 17}


@cache
def polynomial_sections(color, infinity):
    def reverse(rows):
        # A common x power for all components of one section after y=y'/x².
        degree = max(
            max(a.degree(), c.degree() + 2 if c else -1)
            for row in rows
            for a, c in row.values()
        )

        def reflect(poly, n):
            return tuple(int(poly[n - k]) for k in range(n + 1))

        return [
            {
                word: (reflect(a, degree), reflect(c, degree - 2))
                for word, (a, c) in row.items()
            }
            for row in rows
        ]

    families = []
    for rules in presentation(color):
        side = []
        for components in rules:
            records = []
            for h, row in components:
                assert all(not b and not d for a, b, c, d in row.values())
                records.append((h, {w: (a, c) for w, (a, b, c, d) in row.items()}))
            side.append(records)
        families.append(side)
    if not infinity:
        return families
    windows, cuts = families
    converted = [[(h, reverse([r])[0]) for h, r in side] for side in windows]
    converted_cuts = [[] for _ in range(4)]
    for j in range(4):
        for p, r in enumerate(reverse([side[j][1] for side in cuts])):
            converted_cuts[p].append((cuts[p][j][0], r))
    return converted, converted_cuts


def local_model(label, sheet, precision):
    d, a, b, c, e = CENTERS[label]
    ctx.cap = precision
    Jet = gaussian_jets(d, precision)
    radical = Jet(fmpq_poly([]), fmpq_poly([1]))
    u = Jet(fmpq_poly([0, 1]))
    t = Jet.constant(a) + b * radical
    y = sheet * (Jet.constant(c) + e * radical)

    def quartic(t):
        square = t * t
        return -square * square + 14 * square - 1

    if label in ("ram+", "ram-"):
        y = u
        for _ in range(precision.bit_length()):
            t = t - (quartic(t) - u * u) * (-4 * t * t * t + 28 * t).inverse()
    else:
        t = t + u
        half = Jet.constant(2).inverse()
        for _ in range(precision.bit_length()):
            y = (y + quartic(t) * y.inverse()) * half
    assert not (quartic(t) - y * y)

    @cache
    def at(coefficients):
        result = Jet.constant(0)
        for v in reversed(coefficients):
            result = result * t + v
        return result

    def entry(pair):
        a, c = pair
        return at(tuple(map(int, a))) + y * at(tuple(map(int, c)))

    return Jet, entry


def expand(rows, entry):
    return [
        [(h, {w: entry(v) for w, v in row.items()}) for h, row in side] for side in rows
    ]


def certify_disk(label, color, sheet):
    precision = PRECISION[label]
    Jet, entry = local_model(label, sheet, precision)
    windows, cuts = [
        expand(f, entry) for f in polynomial_sections(color, label == "infinity")
    ]
    for p in (0, 1):
        Jet.precision = precision
        losses = []

        def project(n, h, rows, size):
            free, projection, loss = saturated_projection(
                rows, size, DIMENSIONS[n][h], Jet
            )
            losses.append(loss)
            return free, projection

        for _ in quotient_layers(windows, cuts, p, 13, Jet.constant(1), project):
            Jet.precision -= max(losses)
            assert Jet.precision > 0
            losses.clear()
        print(
            label,
            color,
            sheet,
            p,
            "certified; remaining precision",
            Jet.precision,
            flush=True,
        )


def certify_weighted_limit(label, color, sheet, weights):
    precision = 12
    _, entry = local_model(label, sheet, precision)
    windows = expand(polynomial_sections(color, False)[0], entry)
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(101)([-2, 0, 1]))

    def mod(v):
        return int(v.numerator) * pow(int(v.denominator), -1, 101) % 101

    rules = []
    for p, side in enumerate(windows):
        converted = []
        for h, row in side:
            values = {}
            for word, v in row.items():
                weight = sum(weights[(p + k) % 4][a] for k, a in enumerate(word))
                for j in range(precision):
                    a, b = v.real[j], v.imag[j] if v.imag else fmpq(0)
                    if a or b:
                        values[word] = weight + 2 * j, Scalar(field([mod(a), mod(b)]))
                        break
                else:
                    values[word] = weight + 2 * precision, None
            minimum = min(w for w, c in values.values() if c is not None)
            assert all(c is not None or w > minimum for w, c in values.values())
            converted.append(
                (h, {word: c for word, (w, c) in values.items() if w == minimum and c})
            )
        rules.append(converted)
    assert verify_central_profile(rules, Scalar(field.one()))
    print(label, color, sheet, "weighted regular reduction", flush=True)


def large_boundary_factors():
    for q in ((1, 0, 10, 0, 1), (1, 0, 20, 0, -26, 0, 20, 0, 1)):
        for g, e in nmod_poly(list(q), 101).factor()[1]:
            assert e == 1
            yield tuple(map(int, g))


def verify_large_factor(q, color):
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(101)(list(q)))
    t = field.gen() if len(q) > 2 else field(-q[0])

    def quartic(t):
        return -(t**4) + 14 * t * t - 1

    if not quartic(t).is_square():
        field = fq_default_ctx(101, 2 * (len(q) - 1), "a")
        t = fq_default_poly_ctx(field)(list(q)).roots()[0][0]
    y = quartic(t).sqrt()
    for root in (y, -y):

        def entry(pair):
            a, c = pair
            return Scalar(evaluate(a, t, field) + root * evaluate(c, t, field))

        wins, cuts = [expand(f, entry) for f in polynomial_sections(color, False)]

        # At these two factors the original quartet already works.
        def project(n, h, rows, size):
            piv = echelon(rows)
            free = [k for k in range(size) if k not in piv]
            pos = {j: k for k, j in enumerate(free)}
            assert len(free) == DIMENSIONS[n][h]
            return free, [
                {pos[j]: Scalar(field.one())}
                if j in pos
                else {pos[k]: -c for k, c in piv[j].items() if j != k}
                for j in range(size)
            ]

        for start in (0, 1):
            for _ in quotient_layers(
                wins, cuts, start, 13, Scalar(field.one()), project
            ):
                pass
    print("Projective boundary factor", q, color, "certified", flush=True)


WEIGHTS = {
    ("pole+", 0, -1): [[0, 0, -1, -5], [0, -4, -1, -1], [0, 0, -1, -5], [0, 4, 3, 3]],
    ("pole+", 0, 1): [[0, 0, 1, -3], [0, -4, 1, 1], [0, 0, -3, -7], [0, 4, 1, 1]],
    ("pole+", 1, -1): [[0, 0, 1, -3], [0, -4, 1, 1], [0, 0, 1, -3], [0, 4, 5, 5]],
    ("pole+", 1, 1): [[0, 0, 3, -1], [0, -4, 3, 3], [0, 0, -1, -5], [0, 4, 3, 3]],
    ("pole-", 0, 1): [[0, 0, -1, 3], [0, 4, 3, 3], [0, 0, -1, 3], [0, -4, -1, -1]],
    ("pole-", 0, -1): [[0, 0, -3, 1], [0, 4, 1, 1], [0, 0, 1, 5], [0, -4, 1, 1]],
    ("pole-", 1, -1): [[0, 0, -5, -1], [0, 4, -1, -1], [0, 0, -1, 3], [0, -4, -1, -1]],
    ("pole-", 1, 1): [[0, 0, -3, 1], [0, 4, 1, 1], [0, 0, -3, 1], [0, -4, -3, -3]],
}

if __name__ == "__main__":
    for label, color in product(PRECISION, range(2)):
        for sheet in (1, -1) if label == "i" else (1,):
            certify_disk(label, color, sheet)
    for (label, color, sheet), weights in WEIGHTS.items():
        certify_weighted_limit(label, color, sheet, weights)
    for q, color in product(large_boundary_factors(), range(2)):
        verify_large_factor(q, color)
