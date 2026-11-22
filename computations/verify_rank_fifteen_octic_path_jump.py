"""Exact degree jumps of the quadratic Grassmann limits at trace ±1."""

from concurrent.futures import ProcessPoolExecutor
from flint import fmpq_poly
from verify_rank_fifteen_octic_boundary import (
    CYCLOTOMIC,
    IMAGINARY,
    PRECISION,
    SQRT3,
    ZERO,
    grassmann_chart,
    leading_blocks,
    taylor_coefficients,
)
from verify_klein_center_lift import field_paths
from itertools import product


class Cyclotomic:
    """Elements of Q(ζ₂₄) for sparse row reduction."""

    __slots__ = ("value",)

    def __init__(self, value):
        self.value = (
            value % CYCLOTOMIC if isinstance(value, fmpq_poly) else fmpq_poly([value])
        )

    def __bool__(self):
        return bool(self.value)

    def __neg__(self):
        return Cyclotomic(-self.value)

    def __add__(self, other):
        return Cyclotomic(
            self.value + (other.value if isinstance(other, Cyclotomic) else other)
        )

    __radd__ = __add__

    def __mul__(self, other):
        return Cyclotomic(
            self.value * (other.value if isinstance(other, Cyclotomic) else other)
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other.value if isinstance(other, Cyclotomic) else fmpq_poly([other])
        return Cyclotomic(self.value * inverse(other))


def inverse(a):
    gcd, result, _ = a.xgcd(CYCLOTOMIC)
    assert gcd == 1
    return result


def at_ramification(trace, sign_a, sign_b):
    u = sign_a * (SQRT3 if trace == 1 else IMAGINARY)
    v = sign_b * (IMAGINARY if trace == 1 else SQRT3)
    alpha = (v - 1) * inverse(u - 1) % CYCLOTOMIC
    assert (alpha**8 + 14 * alpha**4 + 1) % CYCLOTOMIC == 0
    assert (
        (alpha * alpha + 1) ** 2 + trace * 2 * alpha * (alpha * alpha - 1)
    ) % CYCLOTOMIC == 0
    deriv_inv = inverse((8 * alpha**7 + 56 * alpha**3) % CYCLOTOMIC)
    coeff = taylor_coefficients(alpha)
    scalars = (1, IMAGINARY, SQRT3, SQRT3 * IMAGINARY)

    def leading(coords):
        for n in range(PRECISION):
            letter = n % 2
            order = n // 2
            indices = (2, 3, 6, 7) if letter else (0, 1, 4, 5)
            value = (
                sum(
                    (s * coeff(coords[i], order) for s, i in zip(scalars, indices)),
                    ZERO,
                )
                % CYCLOTOMIC
            )
            if value:
                return n, value * deriv_inv**order % CYCLOTOMIC
        return PRECISION, None

    return leading


def verify_path_jump(job):
    k, a, b, color = job
    rules = []
    for side in leading_blocks(color, at_ramification(k, a, b)):
        rowset = []
        for h, group in side:
            for row in grassmann_chart(
                {cols: Cyclotomic(v) for cols, v in group.items()}, Cyclotomic(1)
            ):
                rowset.append((h, {(x, x ^ h): v for x, v in row.items()}))
        rules.append(rowset)
    basis, _ = field_paths(rules, Cyclotomic(1))
    dimensions = [
        [sum(len(basis[p, n, h]) for h in range(4)) for p in range(4)]
        for n in (3, 4, 5)
    ]
    print(job, dimensions, flush=True)
    assert dimensions == [[20] * 4, [38] * 4, [72] * 4]


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=2) as pool:
        list(pool.map(verify_path_jump, product((1, -1), (-1, 1), (-1, 1), range(2))))
