"""All geometric unit fibers of the four-central-section quotient through degree seven."""

from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx, nmod_poly
from verify_rank_fifteen_bivariate_sixth_cut import (
    B,
    X,
    exceptional_product,
    orbit_coefficients,
    quotient_profile,
    sections,
)
from verify_rank_fifteen_bivariate_word_paths import (
    adjust,
    centers0,
    echelon,
    factors,
    field,
    paths,
    record,
    tensor,
)
from verify_rank_fifteen_quadratic_verticals import split_unit_coordinates
from verify_rank_fifteen_unit_residues import verify_unit_point

fraction = type(field.one)
unit = field.one.n


def primitive(row):
    if not row:
        return {}, unit
    values = iter(row.values())
    content = next(values)
    for value in values:
        if content.total_degree() == 0:
            break
        content = content.gcd(value)
    if content.total_degree() == 0:
        return row, unit
    return {j: z // content for j, z in row.items()}, content


def clear_denominators(row):
    denominator = unit
    for z in row.values():
        denominator *= z.d // denominator.gcd(z.d)
    reduced, content = primitive(
        {j: z.n * (denominator // z.d) for j, z in row.items()}
    )
    return reduced, fraction(content, denominator)


def eliminate(row, scale, pivot, column):
    common = pivot[column].gcd(row[column])
    left, right = pivot[column] // common, row[column] // common
    remainder = {j: z * left for j, z in row.items() if j != column}
    for j, z in pivot.items():
        if j != column:
            value = remainder.get(j, 0) - right * z
            if value:
                remainder[j] = value
            else:
                remainder.pop(j, None)
    reduced, content = primitive(remainder)
    return reduced, scale * fraction(content, left)


def pivot_degree(value, scale, rational_order):
    degree = value.total_degree()
    if rational_order:
        # scale = a/b is reduced. Exactly gcd(value,b) cancels in a*value/b.
        degree += (
            scale.n.total_degree()
            + scale.d.total_degree()
            - 2 * value.gcd(scale.d).total_degree()
        )
    return degree


def primitive_echelon(rows, *, rational_order=False):
    # scale times the primitive polynomial row is the actual rational row.
    pending = [clear_denominators(r) for r in rows if r]
    pivots = []
    while pending:
        width = min(len(row) for row, _ in pending)
        _, j, i = min(
            (pivot_degree(value, scale, rational_order), column, index)
            for index, (row, scale) in enumerate(pending)
            if len(row) == width
            for column, value in row.items()
        )
        pivot, scale = pending.pop(i)
        record((scale * fraction(pivot[j])).n)
        pivots.append((j, pivot))
        for i, (row, scale) in enumerate(pending):
            if j in row:
                pending[i] = eliminate(row, scale, pivot, j)
        pending = [(r, s) for r, s in pending if r]
    result = {}
    for j, row in reversed(pivots):
        for k in list(row):
            if k in result:
                row, _ = eliminate(row, field.one, result[k], k)
        result[j] = row
    return {
        j: {k: fraction(z, row[j]) for k, z in row.items()} for j, row in result.items()
    }


def start_divisors(task):
    start, alternate = task
    # The parent has already certified the universal uncut path and center basis.
    reduce_rows = primitive_echelon if alternate else echelon
    for n, dims in quotient_profile(
        cuts, start, 7, field.one, tensor, adjust, reduce_rows
    ):
        assert sum(dims) == (0, 4, 10, 20, 31, 40, 44, 40)[n]
    assert dims == [10, 10, 10, 10]
    print(
        "Start",
        start,
        "primitive pivots",
        alternate,
        "seventh dimensions",
        dims,
        flush=True,
    )
    return set(factors)


def seventh_divisors():
    global cuts
    basis, image = paths()
    cuts = sections(basis, centers0(basis, image), (0, 1, 3, 4))
    image.cache_clear()
    tasks = [(j, False) for j in range(4)] + [(j, True) for j in (1, 3)]
    # Each fresh worker inherits only the common uncut and center factors.
    with get_context("fork").Pool(processes=4, maxtasksperchild=1) as pool:
        ordinary = pool.map(start_divisors, tasks, chunksize=1)
    return set.union(*ordinary[:4]), set.union(*(ordinary[j] for j in (0, 2, 4, 5)))


def fiber(coefficients, modulus):
    residue = fq_default_ctx(modulus=fmpz_mod_poly_ctx(3)(list(modulus)), var="u")
    x = residue.gen()
    ring = fq_default_poly_ctx(residue)
    common = ring.zero()
    for p in coefficients:
        values = [residue.zero() for _ in range(p.degrees()[1] + 1)]
        for (i, j), z in p.to_dict().items():
            values[j] += int(z) * x**i
        common = common.gcd(ring(values))
    roots = common.roots()
    assert sum(e for _, e in roots) == common.degree()
    return {tuple(map(int, b.to_list())) for b, _ in roots}


def finite_certificate(pair, *, signed_root=False):
    modulus, b_coordinates = pair
    degree = len(modulus) - 1
    for multiple in (1, 2, 3, 4, 6, 12):
        extension = fq_default_ctx(3, degree * multiple, "a")
        ring = fq_default_poly_ctx(extension)
        parameter = ring(list(modulus)).roots()[0][0]
        x = parameter**2 if signed_root else parameter
        if not x.is_square():
            continue
        b = sum(
            (z * parameter**i for i, z in enumerate(b_coordinates)), extension.zero()
        )
        y = x**3 * b
        c = parameter if signed_root else x.sqrt()
        coordinates = split_unit_coordinates(ring.gen(), x, y, c)
        if coordinates is not None:
            verify_unit_point(extension, x, y, coordinates, ((0, 1, 0), (0, 0, 1)))
            print(
                "Certified",
                pair,
                "over extension degree",
                degree * multiple,
                flush=True,
            )
            return
    raise AssertionError(pair)


def horizontal_bound(coefficients, *, signed_root=False):
    # Removing units before elimination cuts down degrees and multiplicities.
    horizontal = X**2 if signed_root else X
    units = (X, X - 1, B, horizontal**3 * B - (1 - horizontal) ** 2)
    if signed_root:
        units += (X + 1,)
    reduced = []
    for p in coefficients:
        for unit_factor in units:
            while not p % unit_factor:
                p = p // unit_factor
        reduced.append(p)
    t = nmod_poly([0, 1], 3)
    bound = None
    for i in (0, 2, 4):
        r = reduced[i].resultant(reduced[i + 1], 1)
        assert r
        univariate = nmod_poly([], 3)
        for (j, k), z in r.to_dict().items():
            assert k == 0
            univariate[j] = int(z)
        for f in (t, t - 1, t + 1) if signed_root else (t, t - 1):
            while not univariate % f:
                univariate = univariate // f
        bound = univariate if bound is None else bound.gcd(univariate)
    return bound


def descent(ordinary, alternate):
    p = exceptional_product(ordinary)
    assert tuple(p.degrees()) == (328, 280)
    curve = X**3 * B + 2 * X + 1
    coefficients = orbit_coefficients(p)
    common = coefficients[0]
    for r in coefficients[1:]:
        common = common.gcd(r)
    assert common == curve
    residual = [r // common for r in coefficients]
    bound = horizontal_bound(residual)
    assert bound.degree() == 70
    pairs = set()
    for f, _ in bound.factor()[1]:
        if f.degree() > 2:
            modulus = tuple(map(int, f))
            pairs.update((modulus, b) for b in fiber(residual, modulus))
    assert len(pairs) == 11

    t = nmod_poly([0, 1], 3)
    q = exceptional_product(alternate)
    assert tuple(q.degrees()) == (391, 347)
    line_bound = t * 0
    for r in orbit_coefficients(q):
        m = r.degrees()[1]
        value = sum(
            (
                int(z) * t ** (a + 3 * (m - b)) * (t - 1) ** b
                for (a, b), z in r.to_dict().items()
            ),
            t * 0,
        )
        line_bound = line_bound.gcd(value)
    assert (
        line_bound
        == (t + 1)
        * (t - 1) ** 39
        * t**175
        * (t * t + t + 2)
        * (t**3 + 2 * t + 2)
        * (t**3 + t * t + t + 2) ** 2
    )
    for f, _ in line_bound.factor()[1]:
        if f.degree() > 2:
            modulus = tuple(map(int, f))
            extension = fq_default_ctx(
                modulus=fmpz_mod_poly_ctx(3)(list(modulus)), var="u"
            )
            x = extension.gen()
            pairs.add((modulus, tuple(map(int, ((x - 1) / x**3).to_list()))))
    assert len(pairs) == 13
    assert sum(len(modulus) - 1 for modulus, _ in pairs) == 44
    return sorted(pairs)


if __name__ == "__main__":
    pairs = descent(*seventh_divisors())
    with ProcessPoolExecutor(max_workers=3) as pool:
        list(pool.map(finite_certificate, pairs))
