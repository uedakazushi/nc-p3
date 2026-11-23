"""Diagonal normalization of the Walsh tensor over a genus-two curve."""

from functools import cache
from itertools import product
from math import prod

import sympy as sp

from flint import fq_default_ctx, fq_default_poly_ctx, nmod_poly

from verify_klein_center_lift import add_multiple, central_equations, field_paths
from verify_rank_fourteen_central_cut import quotient_layers
from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_parameter_orbits import pattern
from verify_rank_fifteen_unit_quartic_towers import encode
from verify_rank_fifteen_unit_residues import verify_unit_point
from verify_rank_fifteen_vertical_all_levels import (
    EllipticUnit,
    ZERO,
    ONE,
    T,
    residue_polynomial,
    sparse_pivot,
)


@cache
def normalize_support(zero_root=None):
    a, b, p, m, q, n, r, l = sp.symbols("a b p m q n r l")
    atoms = (a, b, a + b, p, m, q, n, r, l)
    gamma = pattern(
        a,
        b,
        (p + m) / 2,
        (p - m) / 2,
        (n - q) / 2,
        (n + q) / 2,
        (l - r) / 2,
        (l + r) / 2,
    )
    words, entries = [], []
    for word in product(range(4), repeat=4):
        i, j, k, t = word
        h = i ^ k
        if h != j ^ t:
            continue
        value = sp.expand(
            sum(
                (-1)
                ** (
                    ((i ^ j) & 1) * s
                    + (u >> 1) * (k >> 1)
                    + ((0, 3, 2, 1)[u] >> 1) * (t >> 1)
                )
                * gamma[(h & 2) + s + 4 * u]
                for s, u in product(range(2), range(h & 1, 4, 2))
            )
            / 2
        )
        if value:
            (entry,) = [
                (j, sign)
                for j, atom in enumerate(atoms)
                for sign in (1, -1)
                if value == sign * atom
            ]
            if zero_root is not None and entry[0] == 3 + zero_root:
                continue
            words.append(word)
            entries.append(entry)
    assert len(words) == (36 if zero_root is None else 32)
    incidence = sp.Matrix(
        [[int(w[t] == j) for t, j in product(range(4), repeat=2)] for w in words]
    )
    rows = list(incidence.T.rref()[1])
    cols = list(incidence[rows, :].rref()[1])
    minor = incidence[rows, cols]
    assert len(rows) == 13 and abs(minor.det()) == 4
    coefficients = incidence[:, cols] * minor.inv()
    assert coefficients * incidence[rows, :] == incidence
    assert all(c.is_Integer for c in coefficients)
    normal = {}
    for j, word in enumerate(words):
        atom, sign = entries[j]
        powers = [int(a == atom) for a in range(len(atoms))]
        for index, power in zip(rows, coefficients.row(j), strict=True):
            atom, factor = entries[index]
            powers[atom] -= int(power)
            if power % 2:
                sign *= factor
        assert len(set(p % 2 for p in powers[3:])) == 1
        normal[word] = sign, tuple(powers)
    if zero_root is None:
        assert normal[0, 2, 3, 1][1][3:] == (-1, 1, -1, 1, 1, -1)
    else:
        assert all(all(e % 2 == 0 for e in powers[3:]) for _, powers in normal.values())
    return normal


def normal_windows(base, squares, radical, zero_root=None):
    """Only the six squares and their joint square root enter the coefficients."""
    tensor = {}
    for word, (sign, powers) in normalize_support(zero_root).items():
        odd = powers[3] % 2
        value = sign * prod(x**p for x, p in zip(base, powers[:3], strict=True))
        tensor[word] = (
            value
            * radical**odd
            * prod(
                x ** ((p - odd) // 2) for x, p in zip(squares, powers[3:], strict=True)
            )
        )
    return tensor_windows(tensor)


def tensor_windows(tensor):
    windows = []
    for start in range(4):
        side = []
        for color, k in product(range(4), repeat=2):
            row = {}
            for j in range(4):
                word = (j, j ^ color, k, k ^ color)
                coefficient = tensor.get(word[-start:] + word[:-start])
                if coefficient is not None:
                    row[j, j ^ color] = coefficient
            if row:
                side.append((color, row))
        windows.append(side)
    return windows


class GenusTwo(EllipticUnit):
    radicand = (T**3 - T) ** 2 + 1


def negative_unit_model():
    v = GenusTwo.constant(T)
    i = GenusTwo((ZERO, ONE, ZERO, ZERO))
    root = GenusTwo((ZERO, ZERO, ONE, ZERO))
    return unit_model(v, i, i, root)


def unit_model(v, i, c, root):
    k, h = 1 / v, c * (1 / v - v)
    a, b = type(v).constant(1), (k - 1) / 2
    squares = (
        (k * k - 1 + 2 * h) / 4,
        (k * k - 1 - 2 * h) / 4,
        (k + 1 + h) / 2,
        (k + 1 - h) / 2,
        (k - 1 + h) / 2,
        (k - 1 - h) / 2,
    )
    radical = (-i * c) ** 3 * root * (v**3 - v) ** 2 / v**6
    assert radical**2 == prod(squares)
    return normal_windows((a, b, a + b), squares, radical)


def pivots(rows, unit):
    rows = [r.copy() for r in rows if r]
    echelon = []
    while rows:
        j, k = sparse_pivot(rows)
        row = rows.pop(j)
        inv = row.pop(k).inverse()
        row = {i: v * inv for i, v in row.items()}
        echelon.append((k, row))
        for r in rows:
            if k in r:
                add_multiple(r, row, -r.pop(k))
        rows = [r for r in rows if r]
    out = {}
    for k, row in reversed(echelon):
        for i in list(row):
            if i in out:
                add_multiple(row, out[i], -row[i])
        row[k] = unit
        out[k] = row
    return out


def images(rows, size, unit):
    piv = pivots(rows, unit)
    free = [j for j in range(size) if j not in piv]
    col = {j: k for k, j in enumerate(free)}
    return free, [
        {col[j]: unit}
        if j in col
        else {col[k]: -v for k, v in piv[j].items() if k != j}
        for j in range(size)
    ]


def get_centers(w, unit):
    basis, image = field_paths(w, unit, lambda rows: pivots(rows, unit))
    assert all(
        sum(len(basis[p, n, h]) for h in range(4)) == d
        for p in range(4)
        for n, d in enumerate((1, 4, 10, 20, 35, 56))
    )
    cuts = []
    for h, dim in enumerate((5, 0, 0, 0)):
        rows, size = central_equations(basis, image, h)
        free, proj = images(rows, size, unit)
        assert len(free) == dim
        if h:
            continue
        for j in range(5):
            offset = 0
            parts = []
            for p in range(4):
                words = basis[p, 4, 0]
                parts.append(
                    {
                        word: proj[k + offset][j]
                        for k, word in enumerate(words)
                        if j in proj[k + offset]
                    }
                )
                offset += len(words)
            cuts.append(parts)
    return cuts


def full_fibers(divisor):
    D = T**3 - T
    Y = -(D**2)
    core = divisor // divisor.gcd(Y * (Y - 1))
    gcd = core
    for sign, shift in product((1, -1), range(3)):
        gcd = gcd.gcd(core(sign * T + shift))
    assert gcd.gcd(gcd.derivative()) == ONE
    levels = {}
    for f, e in gcd.factor()[1]:
        assert e == 1
        g = residue_polynomial(Y % f, f)
        levels[g] = levels.get(g, ONE) * f
    for g, f in levels.items():
        poly = nmod_poly(list(g), 3)(Y)
        poly *= pow(int(poly.leading_coefficient()), -1, 3)
        assert f == poly
    return gcd.degree(), tuple(sorted(levels))


def verify_sextic():
    field = fq_default_ctx(3, 12, "a")
    ring = fq_default_poly_ctx(field)
    x = -field.one()
    y = ring([1, 1, 0, 2, 1, 2, 1]).roots()[0][0]
    c = x.sqrt()
    t = ring.gen()
    for v, _ in ((t**3 - t) ** 2 - y / x**3).roots():
        k, h = 1 / v, c * (1 / v - v)
        squares = [k * k - 1 + s * 2 * h for s in (1, -1)]
        squares += [2 * (k + delta + s * h) for delta, s in product((1, -1), repeat=2)]
        if all(z.is_square() for z in squares):
            a, b, p, u, q, r = (z.sqrt() for z in squares)
            coordinates = tuple(map(encode, (k, h, (a + b) / 2, p, u, q, r)))
            verify_unit_point(field, x, y, coordinates, ((0, 1, 0), (0, 0, 1)))
            print("Verified the sextic orbit in extension degree twelve.", flush=True)
            return
    raise AssertionError("No split sextic presentation")


def central_divisor(windows, unit, recorder, degrees):
    centers = get_centers(windows, unit)
    common = recorder.divisor
    cuts = [[(0, centers[j][p]) for j in (1, 2, 3, 4)] for p in range(4)]
    divisor = common
    for start, expected in enumerate(degrees):
        recorder.divisor = common
        if recorder is Recorded:
            Recorded.record.cache_clear()
        for degree in quotient_layers(
            windows,
            cuts,
            start,
            13,
            unit,
            lambda n, h, rows, size: images(rows, size, unit),
        ):
            pass
        assert degree == 13 and recorder.divisor.degree() == expected
        divisor *= recorder.divisor // divisor.gcd(recorder.divisor)
        print(
            "Certified start", start, "with pivot divisor degree", expected, flush=True
        )
    return divisor


def verify():
    Recorded.record.cache_clear()
    Recorded.divisor = ONE
    divisor = central_divisor(
        negative_unit_model(), GenusTwo.constant(1), Recorded, (1049, 2866, 1642, 3847)
    )
    assert divisor.degree() == 8801
    stable, exceptions = full_fibers(divisor)
    assert stable == 156 and len(exceptions) == 8
    assert sorted(len(p) - 1 for p in exceptions) == [2, 2, 2, 3, 3, 4, 4, 6]
    assert [p for p in exceptions if len(p) > 5] == [(1, 1, 0, 2, 1, 2, 1)]
    print(
        "Stable divisor degree",
        stable,
        "and exceptional residue polynomials",
        exceptions,
    )
    verify_sextic()


if __name__ == "__main__":
    verify()
