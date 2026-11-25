"""Two-variable four-central-form cut and descent of the exceptional residue fibers."""

from concurrent.futures import ProcessPoolExecutor
from functools import cache, reduce
from itertools import accumulate, combinations, permutations, product
from operator import xor
from math import prod

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx, nmod_mpoly_ctx, nmod_poly
import sympy as sp

from verify_klein_center_lift import add_multiple
from verify_rank_fifteen_bivariate_word_paths import (
    field,
    K as parameter_k,
    H as parameter_h,
    squares,
    paths,
    centers0,
    adjust,
    word_parity,
    factors,
    tensor,
    echelon,
)


def sections(basis, pivs, selected):
    dims = [len(basis[p, 4, 0]) for p in range(4)]
    offsets = [0, *accumulate(dims)]
    free = [j for j in range(offsets[-1]) if j not in pivs]
    assert len(free) == 5
    out = [[] for p in range(4)]
    for j in selected:
        target = free[j]
        for p in range(4):
            row = {}
            for i, w in enumerate(basis[p, 4, 0]):
                col = offsets[p] + i
                if col == target:
                    row[w] = field.one
                elif col in pivs and target in pivs[col]:
                    row[w] = -pivs[col][target]
            out[p].append(row)
    return out


def quotient_profile(cuts, start, limit, one, entries, power, reduce_rows, *, inspect_degree=None):
    basis = {(0, h): [()] if h == 0 else [] for h in range(4)}
    arrows = {}
    ell = lambda w: word_parity(start, w)

    @cache
    def image(w):
        if not w:
            return {0: one}
        a = w[-1]
        pre = w[:-1]
        h = reduce(xor, pre, 0)
        out = {}
        for j, x in image(pre).items():
            nr = basis[len(pre), h][j]
            e = ell(w) - ell(pre) + ell(nr) - ell(nr + (a,))
            add_multiple(out, arrows[len(pre), h, a][j], x * power(e))
        return out

    for n in range(1, limit + 1):
        for h in range(4):
            words = [w + (a,) for a in range(4) for w in basis[n - 1, h ^ a]]
            off = [sum(w[-1] < a for w in words) for a in range(5)]
            p = (start + n - 2) % 4
            rows = []
            if n >= 2:
                for g, k in product(range(4), repeat=2):
                    for j, pre in enumerate(basis[n - 2, h ^ g]):
                        row = {}
                        shift = None
                        for a in range(4):
                            b = a ^ g
                            tw = (a, b, k, k ^ g)
                            tw = tw[-p:] + tw[:-p]
                            if tw not in entries:
                                continue
                            exponent, entry = entries[tw]
                            for idx, z in arrows[n - 2, h ^ g, a][j].items():
                                nr = basis[n - 1, h ^ b][idx]
                                e = (
                                    exponent
                                    - ell(pre + (a,))
                                    + ell(nr)
                                    - ell(nr + (b,))
                                )
                                if shift is None:
                                    shift = e % 2
                                col = off[b] + idx
                                row[col] = row.get(col, 0) + entry * z * power(
                                    e - shift
                                )
                        rows.append({j: z for j, z in row.items() if z})
            if n == 4 and h == 0:
                for sect in cuts[start]:
                    row = {}
                    for w, x in sect.items():
                        pre = w[:-1]
                        b = w[-1]
                        for j, z in image(pre).items():
                            nr = basis[3, b][j]
                            e = ell(w) - ell(pre) + ell(nr) - ell(nr + (b,))
                            col = off[b] + j
                            row[col] = row.get(col, 0) + x * z * power(e)
                    rows.append({j: z for j, z in row.items() if z})
            piv = reduce_rows(rows)
            free = [j for j in range(len(words)) if j not in piv]
            ids = {j: k for k, j in enumerate(free)}
            basis[n, h] = [words[j] for j in free]
            proj = [
                {ids[j]: one}
                if j in ids
                else {ids[k]: -z for k, z in piv[j].items() if k != j}
                for j in range(len(words))
            ]
            for a in range(4):
                arrows[n - 1, h ^ a, a] = proj[off[a] : off[a + 1]]
        dims = [len(basis[n, h]) for h in range(4)]
        if inspect_degree is not None:
            inspect_degree(n, basis, image)
        yield n, dims
    image.cache_clear()


def quotient(cuts, start, limit):
    for n, dims in quotient_profile(cuts, start, limit, field.one, tensor, adjust, echelon):
        assert sum(dims) == (1, 4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1, 0)[n], (
            n,
            dims,
        )
        print("cut", start, n, dims, len(factors), flush=True)
        yield n, set(factors)


allowed = [parameter_k.n, parameter_h.n, (parameter_k - 1).n, (parameter_k + 1).n] + [
    p for v in squares for p, e in v.n.factor()[1]
]


def strip(f):
    for p in allowed:
        while f and not f % p:
            f = f // p
    return f


def fullrank(rows):
    rs = [r.copy() for r in rows]
    while rs:
        candidates = sorted(
            ((i, j) for i, r in enumerate(rs) for j in r),
            key=lambda ij: (len(rs[ij[0]]), rs[ij[0]][ij[1]].n.total_degree(), ij[1]),
        )
        for i, j in candidates:
            if strip(rs[i][j].n).total_degree() == 0:
                break
        else:
            break
        row = rs.pop(i)
        inv = field.one / row.pop(j)
        row = {k: z * inv for k, z in row.items()}
        for r in rs:
            if j in r:
                add_multiple(r, row, -r.pop(j))
    if not rs:
        return
    cols = sorted(set().union(*(r.keys() for r in rs)))
    n = len(rs)
    # Test the remaining maximal minors in the localization at the chart.
    k, h = sp.symbols("k h")

    def sym(f):
        return sum(int(z) * k ** e[0] * h ** e[1] for e, z in f.to_dict().items())

    signs = [
        (perm, (-1) ** sum(a > b for j, a in enumerate(perm) for b in perm[j + 1 :]))
        for perm in permutations(range(n))
    ]
    gens = []
    for subset in combinations(cols, n):
        det = sum(
            (
                sg * prod(row.get(subset[j], field.zero) for row, j in zip(rs, perm))
                for perm, sg in signs
            ),
            field.zero,
        )
        if not det:
            continue
        f = strip(det.n)
        if f.total_degree() == 0:
            return
        gens.append(sym(f))
        gb = sp.groebner(gens, k, h, modulus=3, order="grevlex")
        if gb.contains(sp.Integer(1)):
            return
        gens = list(gb)
    # The input chart excludes the zeroes removed during row reduction.
    w = sp.symbols("w")
    C0 = prod(map(sym, allowed))
    gb = sp.groebner([*gens, 1 - w * C0], w, k, h, modulus=3, order="grevlex")
    assert gb.contains(sp.Integer(1)), gb


def products(basis, image, cs):
    for p in range(4):
        for arrow in (None, 0, 1, 2, 3):
            rows = []
            for section in cs[p]:
                out = {}
                for w, x in section.items():
                    full = w if arrow is None else w + (arrow,)
                    delta = (
                        0
                        if arrow is None
                        else word_parity(p, w)
                        - word_parity(p, full)
                        - word_parity(p, (arrow,))
                    )
                    add_multiple(out, image(p, full), x * adjust(delta))
                rows.append(out)
            fullrank(rows)


def start_certificate(start):
    factors.clear()
    basis, image = paths()
    pivots = centers0(basis, image)
    if start == 0:
        products(basis, image, sections(basis, pivots, range(5)))
    cuts = sections(basis, pivots, (0, 1, 3, 4))
    image.cache_clear()
    through = {}
    for n, divisors in quotient(cuts, start, 6):
        if n >= 5:
            through[n] = divisors
    return through


cover = nmod_mpoly_ctx.get(("v", "c"), 3)
v, c = cover.gens()
base = nmod_mpoly_ctx.get(("X", "B"), 3)
X, B = base.gens()
K, H = sp.symbols("K H")


def pullback(source):
    p = sp.Poly(source.replace("^", "**"), K, H, modulus=3)
    degree = p.total_degree()
    return sum(
        (int(a) % 3) * v ** (degree - i - j) * (c * (1 - v * v)) ** j
        for (i, j), a in p.terms()
    )


def exceptional_product(factor_set):
    chart = (
        K,
        H,
        K + 1,
        K - 1,
        K + 1 + H,
        K + 1 - H,
        K - 1 + H,
        K - 1 - H,
        K * K - 1 + H,
        K * K - 1 - H,
    )
    excluded = v * c * (c * c - 1)
    for f in chart:
        excluded *= pullback(str(f))
    result = cover.constant(1)
    for source in sorted(factor_set):
        f = pullback(source)
        while (common := f.gcd(excluded)).total_degree() > 0:
            f = f // common
        result *= f // result.gcd(f)
    return result


@cache
def reduced_power(n):
    if n < 6:
        return tuple(base.constant(int(j == n)) for j in range(6))
    return tuple(
        B * g + 2 * t - r
        for g, t, r in zip(
            reduced_power(n - 6),
            reduced_power(n - 2),
            reduced_power(n - 4),
            strict=True,
        )
    )


def orbit_coefficients(polynomial):
    coefficients = [[base.constant(0) for j in range(6)] for i in range(2)]
    for (j, i), a in polynomial.to_dict().items():
        for e, r in enumerate(reduced_power(j)):
            coefficients[i % 2][e] += int(a) * X ** (i // 2) * r
    return [r for row in coefficients for r in row if r]


def horizontal_bound(polynomial):
    nonzero = orbit_coefficients(polynomial)
    # Every bad pair kills every coefficient. Any subset of pairwise resultants
    # therefore gives a necessary horizontal condition, including at a degree drop.
    bound = None
    for i in range(2):
        for j in range(i + 1, len(nonzero)):
            resultant = nonzero[i].resultant(nonzero[j], 1)
            if not resultant:
                continue
            poly = nmod_poly([], 3)
            for (xpower, bpower), a in resultant.to_dict().items():
                assert bpower == 0
                poly[xpower] = int(a)
            bound = poly if bound is None else bound.gcd(poly)
            if bound.degree() == 0 or (i == 1 and j >= min(5, len(nonzero) - 1)):
                return bound * pow(int(bound.leading_coefficient()), -1, 3)
    assert bound is not None
    return bound * pow(int(bound.leading_coefficient()), -1, 3)


def sixth_fiber(polynomial):
    f = fq_default_ctx(modulus=fmpz_mod_poly_ctx(3)([2, 1, 1]), var="t")
    x = f.gen()
    ring = fq_default_poly_ctx(f)
    common = ring.zero()
    for coefficient in orbit_coefficients(polynomial):
        values = [f.zero() for _ in range(coefficient.degrees()[1] + 1)]
        for (i, j), z in coefficient.to_dict().items():
            values[j] += int(z) * x**i
        common = common.gcd(ring(values))
    assert common == ring.gen() - (x + 1)
    print("The simultaneous exceptional fiber is B = X + 1.", flush=True)


def verify_orbits(results):
    t = nmod_poly([0, 1], 3)
    for n, degree, expected in (
        (5, (7, 7), t**0),
        (6, (55, 44), t**136 * (t - 1) ** 9 * (t * t + t + 2)),
    ):
        polynomial = exceptional_product(set.union(*(r[n] for r in results)))
        assert tuple(polynomial.degrees()) == degree
        bound = horizontal_bound(polynomial)
        assert bound == expected
        print("Length", n, "horizontal bound factors", bound.factor()[1], flush=True)
        if n == 6:
            sixth_fiber(polynomial)


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=4) as workers:
        results = list(workers.map(start_certificate, range(4)))
    verify_orbits(results)
