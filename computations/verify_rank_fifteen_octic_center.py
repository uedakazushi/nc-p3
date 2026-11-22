"""Uniform quartic center and residual closed-point certificates on the genus-three octic."""

from functools import cache
from itertools import product

from flint import fmpz_poly, nmod_poly, nmod
from sympy import primerange, sqrt_mod

from rank_fifteen_octic import finite_boundary, octic_frame, finite_windows
from rank_fifteen_octic_field import UNIT, windows
from rank_fourteen_function_field import ONE, irreducible_factors
from verify_klein_center_lift import add_multiple, central_equations, field_paths
from verify_rank_fourteen_geometric_reduction import verify_central_profile


def residue(value):
    # t,Y,i,sqrt(3)=(2,46,33,49) at 109 select independent rows and pivot columns.
    def evaluate(v):
        a, b, c, d = v.coordinates
        return (
            int((a(2) + 33 * b(2) + 46 * (c(2) + 33 * d(2))) % 109)
            * pow(int(v.denominator(2)) % 109, -1, 109)
            % 109
        )

    return (evaluate(value.a) + (49 * evaluate(value.b) if value.b else 0)) % 109


def chosen_rows(rows, central):
    pivots = {}
    order = (
        sorted(range(len(rows)), key=lambda j: (len(rows[j]), j))
        if central
        else range(len(rows))
    )
    for j in order:
        row = {k: v for k, c in rows[j].items() if (v := residue(c))}
        for k, previous in pivots.items():
            if k not in row:
                continue
            factor = row.pop(k)
            for a, b in previous.items():
                value = (row.get(a, 0) - factor * b) % 109
                if value:
                    row[a] = value
                else:
                    row.pop(a, None)
        if row:
            k = min(row)
            inverse = pow(row.pop(k), -1, 109)
            pivots[k] = {a: x * inverse % 109 for a, x in row.items()}
            yield rows[j].copy(), k


def exceptional_support(color):
    rules = windows(color)
    bad = set()

    @cache
    def record(poly):
        bad.update(irreducible_factors(fmpz_poly(list(poly))))

    def store(poly):
        record(tuple(map(int, poly)))

    def denominators(value):
        store(value.a.denominator)
        if value.b:
            store(value.b.denominator)

    for side in rules:
        for _, row in side:
            for value in row.values():
                denominators(value)

    def reduce(rows, central=False):
        pivots = {}
        num = den = ONE
        for row, k in chosen_rows(rows, central):
            for a, previous in pivots.items():
                if a in row:
                    add_multiple(row, previous, -row.pop(a))
            coefficient = row.pop(k)
            if central:
                base = coefficient.norm()
                u, g = base.norm_numerator(), base.denominator**4
                left, right = num.gcd(g), den.gcd(u)
                num = (num // left) * (u // right)
                den = (den // right) * (g // left)
            inverse = coefficient.inverse()
            pivots[k] = {a: v * inverse for a, v in row.items()}
        if central:
            store(num)
            store(den)
            return pivots
        for k in reversed(pivots):
            row = pivots[k]
            for a in list(row):
                if a in pivots:
                    add_multiple(row, pivots[a], -row.pop(a))
        for k, row in pivots.items():
            row[k] = UNIT
            for value in row.values():
                denominators(value)
        return pivots

    basis, image = field_paths(rules, UNIT, reduce)
    assert all(
        sum(len(basis[p, n, h]) for h in range(4)) == d
        for p in range(4)
        for n, d in enumerate((1, 4, 10, 20, 35, 56))
    )
    for h, d in enumerate((3, 2, 0, 0)):
        rows, width = central_equations(basis, image, h)
        assert width - len(reduce(rows, True)) == d
    return sorted(bad - finite_boundary())


def witness(q):
    boundary = finite_boundary()
    for p in primerange(97, 2000):
        if p % 12 != 1:
            continue
        g = nmod_poly(list(q), p)
        if g.degree() != len(q) - 1:
            continue
        for root, e in g.roots():
            t = int(root)
            if e != 1 or not all(nmod_poly(list(f), p)(t) for f in boundary):
                continue
            if len(sqrt_mod(t**8 + 14 * t**4 + 1, p, all_roots=True)) == 2:
                return p, t
    raise ValueError(q)


@cache
def certify(color, p, root):
    t = nmod(root, p)
    js = [j for j in range(p) if (j**4 - j * j + 1) % p == 0]
    ys = sqrt_mod(int(t**8 + 14 * t**4 + 1), p, all_roots=True)
    assert len(js) == 4 and len(ys) == 2
    unit = nmod(1, p)
    for j, y in product(js, ys):
        frame = octic_frame(t, nmod(y, p), nmod(j, p), color)
        zz, ll, *_ = frame
        assert (
            zz
            * ll
            * (zz * zz - 1)
            * (ll * ll - 1)
            * (ll * ll * zz - 1)
            * (ll * ll * zz * zz - 1)
        )
        assert verify_central_profile(finite_windows(frame, unit, lambda x: x), unit)


def run(color):
    support = exceptional_support(color)
    assert len(support) == 17 and sum(len(f) - 1 for f in support) == 888
    for q in support:
        p, t = witness(q)
        certify(color, p, t)
        print(
            "Presentation",
            color,
            "degree",
            len(q) - 1,
            "split witness",
            p,
            t,
            flush=True,
        )


if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(max_workers=2) as pool:
        list(pool.map(run, (0, 1)))
