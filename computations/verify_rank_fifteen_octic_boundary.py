"""Certify the genus-three relation lattices above trace 0, 2 and -2 at 109."""

from concurrent.futures import ProcessPoolExecutor
from functools import cache
from itertools import combinations
from math import comb

from flint import fmpq_poly, fq_default_ctx, fq_default_poly_ctx

from rank_fifteen_octic_field import UNIT, windows
from rank_fourteen_function_field import ONE
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile

CYCLOTOMIC = fmpq_poly([1, 0, 0, 0, -1, 0, 0, 0, 1])
ETA = fmpq_poly([0, 1])
IMAGINARY = ETA**6
SQRT2 = (ETA**3 + ETA**21) % CYCLOTOMIC
SQRT3 = (ETA**2 + ETA**22) % CYCLOTOMIC
ZERO = fmpq_poly([])
OCTIC = (1, 0, 0, 0, 14, 0, 0, 0, 1)
PRECISION = 10


def clear_coordinates(values):
    """Use one denominator for the entire projective coordinate tuple."""
    denominator = ONE
    for value in values.values():
        for part in (value.a, value.b):
            if part:
                denominator *= part.denominator // denominator.gcd(part.denominator)

    def components(part):
        if not part:
            return [()] * 4
        scale = denominator // part.denominator
        return [tuple(map(int, scale * poly)) for poly in part.coordinates]

    return {
        columns: components(value.a) + components(value.b)
        for columns, value in values.items()
    }


@cache
def pluecker_coordinates(color):
    result = []
    zero = UNIT.lift(0)
    for side in windows(color):
        blocks = []
        for character in range(4):
            rows = [row for tag, row in side if tag == character]
            assert len(rows) in (1, 2)

            def entry(row, col):
                return row.get((col, col ^ character), zero)

            values = {}
            for columns in combinations(range(4), len(rows)):
                if len(rows) == 1:
                    value = entry(rows[0], columns[0])
                else:
                    a, b = columns
                    value = entry(rows[0], a) * entry(rows[1], b) - entry(
                        rows[0], b
                    ) * entry(rows[1], a)
                if value:
                    values[columns] = value
            blocks.append((character, clear_coordinates(values)))
        result.append(blocks)
    return result


def taylor_coefficients(alpha):
    @cache
    def alpha_power(n):
        return alpha**n % CYCLOTOMIC

    @cache
    def coefficient(poly, n):
        return (
            sum(
                (
                    c * comb(j, n) * alpha_power(j - n)
                    for j, c in enumerate(poly)
                    if j >= n and c
                ),
                ZERO,
            )
            % CYCLOTOMIC
        )

    return coefficient


def taylor_evaluation(trace, t_sign):
    alpha = t_sign * IMAGINARY if trace == 0 else -(trace // 2) + t_sign * SQRT2
    height = IMAGINARY if trace == 0 else SQRT3
    y0 = 2 * alpha * (alpha**2 - 1) * height % CYCLOTOMIC

    coefficient = taylor_coefficients(alpha)
    assert y0**2 % CYCLOTOMIC == coefficient(OCTIC, 0)
    gcd, inverse, _ = (2 * y0).xgcd(CYCLOTOMIC)
    assert gcd == 1
    assert all(int(c.denominator) % 109 for c in inverse)
    y = [y0]
    for n in range(1, PRECISION):
        square = sum((y[k] * y[n - k] for k in range(1, n)), ZERO)
        y.append((coefficient(OCTIC, n) - square) * inverse % CYCLOTOMIC)

    scalars = (1, IMAGINARY, SQRT3, IMAGINARY * SQRT3)

    def leading(components):
        def part(indices, n):
            return (
                sum(
                    (
                        s * coefficient(components[i], n)
                        for s, i in zip(scalars, indices)
                    ),
                    ZERO,
                )
                % CYCLOTOMIC
            )

        ordinate = []
        for n in range(PRECISION):
            ordinate.append(part((2, 3, 6, 7), n))
            value = (
                part((0, 1, 4, 5), n)
                + sum((y[n - j] * ordinate[j] for j in range(n + 1)), ZERO)
            ) % CYCLOTOMIC
            if value:
                return n, value
        return PRECISION, None

    return leading


def leading_blocks(color, leading):
    result = []
    for side in pluecker_coordinates(color):
        blocks = []
        for character, coordinates in side:
            values = {columns: leading(parts) for columns, parts in coordinates.items()}
            order = min(n for n, _ in values.values())
            assert order < PRECISION
            blocks.append(
                (character, {cols: v for cols, (n, v) in values.items() if n == order})
            )
        result.append(blocks)
    return result


def grassmann_chart(group, unit):
    pivot = next(iter(group))
    for j, col in enumerate(pivot):
        row = {col: unit}
        for k in range(4):
            if k in pivot:
                continue
            cols = list(pivot)
            cols[j] = k
            key = tuple(sorted(cols))
            if key in group:
                inversions = sum(
                    a > b for i, a in enumerate(cols) for b in cols[i + 1 :]
                )
                row[k] = (-1) ** inversions * group[key] / group[pivot]
        yield row


def reduce_blocks(blocks, field, root):
    def evaluate(poly):
        value = field.zero()
        for c in reversed(list(poly)):
            n, d = int(c.numerator), int(c.denominator)
            assert d % 109
            value = value * root + n * pow(d, -1, 109)
        return value

    result = []
    for side in blocks:
        rules = []
        for character, coords in side:
            group = {
                cols: value
                for cols, poly in coords.items()
                if not (value := evaluate(poly)).is_zero()
            }
            assert group
            for row in grassmann_chart(group, field.one()):
                rules.append(
                    (character, {(a, a ^ character): Scalar(c) for a, c in row.items()})
                )
        result.append(rules)
    return result


def verify(job):
    trace, sign, color = job
    blocks = leading_blocks(color, taylor_evaluation(trace, sign))
    field = fq_default_ctx(109, 2, "a")
    roots = [
        v for v, _ in fq_default_poly_ctx(field)(list(map(int, CYCLOTOMIC))).roots()
    ]
    assert len(roots) == 8
    for j in (8, 41):
        root = next(v for v in roots if v * v == j)
        assert verify_central_profile(
            reduce_blocks(blocks, field, root), Scalar(field.one())
        )
    print("Verified boundary", trace, sign, color, flush=True)


if __name__ == "__main__":
    assert IMAGINARY**2 % CYCLOTOMIC == -1
    assert SQRT2**2 % CYCLOTOMIC == 2 and SQRT3**2 % CYCLOTOMIC == 3
    assert (SQRT3 + IMAGINARY - 2 * ETA**2) % CYCLOTOMIC == 0
    jobs = [
        (k, s, h)
        for k in (0, 2, -2)
        for s in ((-1, 1) if k == 0 else (1,))
        for h in range(2)
    ]
    with ProcessPoolExecutor(max_workers=2) as pool:
        list(pool.map(verify, jobs))
