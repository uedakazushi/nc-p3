"""All prime-field points with y nonzero on the regularity plane at 3, 5 and 7."""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from itertools import product

from flint import fq_default_ctx, fq_default_poly_ctx

from verify_rank_fifteen_parameter_orbits import amplitudes
from verify_rank_fifteen_pencil_linear_residues import graded_profile
from verify_rank_fifteen_regular_orbits import square_roots, tensor_windows
from verify_rank_fourteen_geometric_reduction import Scalar, verify_central_profile


@cache
def contexts(prime, degree):
    field = fq_default_ctx(prime, degree, "a")
    return field, fq_default_poly_ctx(field)


def lift_coordinates(field, polynomials, x, y):
    for v, _ in polynomials([-y, x * x, -8 * x, 16]).roots():
        U, V = x / v - 3, y / (v * v)
        if not (U.is_square() and V.is_square()):
            continue
        for K, H in product(square_roots(U), square_roots(V)):
            dr, dw = (K + H) ** 2 - 1, (K - H) ** 2 - 1
            if not (dr.is_square() and dw.is_square()):
                continue
            rs = [K + H + h for h in square_roots(dr)]
            ws = [K - H + h for h in square_roots(dw)]
            ns = [
                u.sqrt() for u, _ in polynomials([V, 1 - U, 1]).roots() if u.is_square()
            ]
            for r, w, n in product(rs, ws, ns):
                if r.is_square() and w.is_square():
                    s, t = r.sqrt(), w.sqrt()
                    assert (r + 1 / r + w + 1 / w) / 4 == K and (
                        r + 1 / r - w - 1 / w
                    ) / 4 == H
                    assert V * (U + 3) / (U - 1) ** 2 == x and V**3 / (U - 1) ** 4 == y
                    assert n**4 + (1 - U) * n * n + V == field.zero()
                    return s, t, n
    return None


def verify_point(args):
    prime, x, y = args
    for degree in (1, 2, 3, 4, 6, 8, 12, 16, 24):
        field, polynomials = contexts(prime, degree)
        coords = lift_coordinates(field, polynomials, field(x), field(y))
        if coords is None:
            continue
        s, t, n = coords
        windows = tensor_windows(amplitudes(s, t / s, n, field.one()), field)
        unit = Scalar(field.one())
        pairs = (
            ((1, 0, 1), (0, 1, 1))
            if (prime, x, y) == (7, 0, 2)
            else ((1, 3, 5), (2, 7, 11))
        )
        try:
            assert graded_profile(windows, unit, pairs)
        except AssertionError:
            assert verify_central_profile(windows, unit)
        return degree
    raise AssertionError((prime, x, y))


def verify():
    with ProcessPoolExecutor(max_workers=3) as pool:
        for prime in (3, 5, 7):
            points = [(prime, x, y) for x in range(prime) for y in range(1, prime)]
            counts = Counter(pool.map(verify_point, points))
            assert sum(counts.values()) == prime * (prime - 1)
            print(
                prime,
                len(points),
                "quotient residues certified; witness degrees",
                sorted(counts.items()),
                flush=True,
            )


if __name__ == "__main__":
    verify()
