"""The two prime-field pencils Y/(1-X)^2 on the ordinary unit chart."""

from flint import fq_default_ctx, fq_default_poly_ctx
from verify_rank_fifteen_genus_two_unit_line import (
    GenusTwo,
    ONE,
    ZERO,
    T,
    Recorded,
    central_divisor,
    full_fibers,
    normal_windows,
    unit_model,
)
from verify_rank_fifteen_quadratic_verticals import split_unit_coordinates
from verify_rank_fifteen_unit_residues import verify_unit_point


def branch_windows():
    c = GenusTwo.constant(T)
    k = c
    h = k * k - 1
    unit = GenusTwo.constant(1)
    b = (k - 1) / 2
    squares = (
        (k * k - 1 + 2 * h) / 4,
        (k * k - 1 - 2 * h) / 4,
        (k + 1 + h) / 2,
        (k + 1 - h) / 2,
        (k - 1 + h) / 2,
        (k - 1 - h) / 2,
    )
    assert not squares[0] and all(squares[1:])
    return normal_windows((unit, b, unit + b), squares, unit, zero_root=0)


def negative_windows():
    q = GenusTwo.constant(T)
    i = GenusTwo((ZERO, ONE, ZERO, ZERO))
    a = -i
    b = (q - q**3) / (i - a)
    v = q + a * b
    d = b**3 - b
    assert v**3 - v == i * d and b * b == -((q**3 - q) ** 2)
    return unit_model(v, i, 1 / b, d)


def exceptional(g):
    d = len(g) - 1
    if d <= 2:
        return
    field = fq_default_ctx(3, 8 if d == 4 else 12, "a")
    ring = fq_default_poly_ctx(field)
    x = 1 / ring(list(g)).roots()[0][0]
    y = -((1 - x) ** 2)
    coordinates = split_unit_coordinates(ring.gen(), x, y, x.sqrt())
    assert coordinates is not None
    verify_unit_point(field, x, y, coordinates, ((0, 1, 0), (0, 0, 1)))
    print("Verified the reciprocal horizontal residue polynomial", g, flush=True)


def verify():
    for branch in (True, False):
        Recorded.record.cache_clear()
        Recorded.divisor = ONE
        windows = branch_windows() if branch else negative_windows()
        degrees = (7, 3, 5, 5) if branch else (1174, 1285, 1375, 1074)
        divisor = central_divisor(windows, GenusTwo.constant(1), Recorded, degrees)
        if branch:
            expected = T * (T - 1) * (T + 1) * (T * T + 1)
            assert divisor == expected * (T * T + T + 2)
            assert divisor.gcd(divisor(-T)) == expected
        else:
            assert divisor.degree() == 4438
            stable, levels = full_fibers(divisor)
            assert stable == 126 and sorted(len(g) - 1 for g in levels) == [
                1,
                2,
                2,
                2,
                3,
                3,
                4,
                4,
            ]
            for g in levels:
                exceptional(g)


if __name__ == "__main__":
    verify()
