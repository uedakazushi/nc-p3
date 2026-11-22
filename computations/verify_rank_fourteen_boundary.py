"""Classify the projective tensor limits in the homogeneous rank-fourteen frame."""

from itertools import product

from flint import fmpz_mod_poly_ctx, fmpz_poly, fq_default_ctx, fq_default_poly_ctx
from rank_fourteen_function_field import (
    imaginary,
    irreducible_factors,
    one,
    ordinate,
    parameter,
)
from sympy import QQ_I
from verify_klein_center_lift import echelon
from verify_rank_fourteen_cyclotomic import boundary
from verify_rank_fourteen_descent import frame, slopes
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile


def homogeneous_tensor():
    t = parameter
    r, lam = slopes(t)
    z = r * r
    x, b, e, f = frame(ordinate / ((t + 1) * (t + 3)), imaginary)
    X = ((one, one / lam), (one / lam, z))
    Y = ((one, lam), (lam * z, one))

    def multiply(left, right):
        return [
            [sum(left[j][k] * right[k][l] for k in range(2)) for l in range(2)]
            for j in range(2)
        ]

    def thick(matrix, scale):
        (a, b), (c, d) = matrix
        det = a * d - b * c
        inverse = ((d / det, -b / det), (-c / det, a / det))
        lower = multiply(multiply(Y, inverse), X)
        return [[*matrix[j], *(scale * v for v in X[j])] for j in range(2)] + [
            [*(v / scale for v in Y[j]), *lower[j]] for j in range(2)
        ]

    blocks = (
        thick(((r * x, b), (z * b, r * x)), one),
        list(zip(*thick(((r * e, f), (z * f, r * e)), lam))),
    )

    def thin(a, b):
        return z if a == 3 and b == 1 else one

    tensor = {}
    for a, b, c in product(range(4), repeat=3):
        d = a ^ b ^ c
        tensor[a, b, c, d] = (
            blocks[a ^ b][a][c] if a ^ b < 2 else thin(a, b) * thin(c, d ^ 1)
        )
    return tensor


def flattenings(tensor):
    rules = []
    for p in range(4):
        rows = [{} for _ in range(16)]
        for indices, value in tensor.items():
            a, b, c, d = indices[p:] + indices[:p]
            if value:
                rows[4 * c + d][4 * a + b] = value
        side = []
        for row in echelon(rows).values():
            entries = {divmod(k, 4): v for k, v in row.items()}
            a, b = next(iter(entries))
            assert all((i ^ j) == (a ^ b) for i, j in entries)
            side.append((a ^ b, entries))
        rules.append(side)
    rows = [{} for _ in range(16)]
    for (a, b, c, d), value in tensor.items():
        if value:
            rows[4 * b + d][4 * a + c] = value
    return rules, len(echelon(rows))


def specialize(tensor, t, i, y, wrap=lambda x: x):
    def poly(value):
        result = 0 * t
        for c in reversed(list(value)):
            result = result * t + int(c)
        return result

    def entry(value):
        a, b, c, d = map(poly, value.coordinates)
        return wrap((a + i * b + y * (c + i * d)) / poly(value.denominator))

    return {word: entry(c) for word, c in tensor.items()}


def test_regular(tensor, f):
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(101)(list(f)), var="a")
    t = field.gen()
    cubic = 2 * (t - 1) * (t + 1) * (t + 3)
    if not cubic.is_square():
        field = fq_default_ctx(101, 2 * (len(f) - 1), "a")
        t = fq_default_poly_ctx(field)(list(f)).roots()[0][0]
        cubic = 2 * (t - 1) * (t + 1) * (t + 3)
    y = cubic.sqrt()
    for i in (10, 91):
        for root in (y, -y):
            values = specialize(tensor, t, field(i), root, Scalar)
            rules, opposite = flattenings(values)
            assert list(map(len, rules)) == [6] * 4 and opposite == 14
            for p in range(4):
                assert len({ix[p] for ix, z in values.items() if z}) == 4
            assert verify_central_profile(rules, Scalar(field.one()))
    print("Regular boundary residues", f, flush=True)


def leading(poly, t0):
    if not poly:
        return None
    if t0 is None:
        return -poly.degree(), int(poly.leading_coefficient())
    order = 0
    factor = fmpz_poly([-t0, 1])
    while not poly(t0):
        poly //= factor
        order += 1
    return order, int(poly(t0))


def ramified_limit(tensor, t0):
    coefficients = {}
    for word, value in tensor.items():
        vq, lc = leading(value.denominator, t0)
        pairs = [leading(g, t0) for g in value.coordinates]
        terms = []
        for sheet in range(2):
            pair = pairs[2 * sheet : 2 * sheet + 2]
            if any(term is not None for term in pair):
                v = min(term[0] for term in pair if term is not None)
                a, b = (
                    term[1] if term is not None and term[0] == v else 0 for term in pair
                )
                order = 2 * (v - vq) + sheet * (1 if t0 is not None else -3)
                terms.append((order, QQ_I(a, b) / QQ_I(lc)))
        if terms:
            coefficients[word] = min(terms, key=lambda term: term[0])
    minimum = min(order for order, c in coefficients.values())
    return {word: c for word, (order, c) in coefficients.items() if order == minimum}


def cubic_dimensions(rules):
    dims = []
    for p in range(4):
        rows = []
        for a in range(4):
            rows.extend(
                [
                    {16 * i + 4 * j + a: c for (i, j), c in row.items()}
                    for h, row in rules[p]
                ]
            )
            rows.extend(
                [
                    {16 * a + 4 * i + j: c for (i, j), c in row.items()}
                    for h, row in rules[(p + 1) % 4]
                ]
            )
        dims.append(64 - len(echelon(rows)))
    return dims


if __name__ == "__main__":
    tensor = homogeneous_tensor()
    assert set().union(
        *(irreducible_factors(v.denominator) for v in tensor.values())
    ) == {(-1, 1), (1, 1), (3, 1)}
    special = {(1, 1), (3, 1), (100, 1), (21, 1), (82, 1)}
    regular = boundary(101) - special
    assert regular == {(94, 99, 1), (1, 6, 1), (40, 81, 1), (59, 20, 1)}
    for f in sorted(regular):
        test_regular(tensor, f)
    i = QQ_I(0, 1)
    t = -1 + 2 * i
    for component in (i, -i):
        for s in (1, -1):
            rules, opp = flattenings(
                specialize(tensor, t, component, (t + 1) * (t + 3) * s)
            )
            assert list(map(len, rules)) == [4] * 4 and opp == 4
    for t0, adjacent, expected in (
        (-3, [2, 1, 2, 1], 2),
        (-1, [2] * 4, 1),
        (1, [6] * 4, 6),
        (None, [2] * 4, 1),
    ):
        rules, opp = flattenings(ramified_limit(tensor, t0))
        assert list(map(len, rules)) == adjacent and opp == expected
        if t0 == 1:
            assert cubic_dimensions(rules) == [28, 28, 22, 22]
    print(
        "All 48 geometric boundary parameters classified in the projective tensor frame.",
        flush=True,
    )
