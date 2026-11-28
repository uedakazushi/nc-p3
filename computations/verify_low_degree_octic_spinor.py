"""Verify the numerical alternatives in the spinor obstruction to low-degree factors."""
from sympy import expand, symbols


def candidates(degree, bound, ramified_vertex=False, plane=False):
    result = []
    for n in range(bound + 1):
        h_square = 4 if ramified_vertex else degree
        length = 2 * n + 1 if ramified_vertex else n
        spinor_marks = 2 * n if ramified_vertex else n
        c2 = 6 * h_square - length
        intersection = 2 * length - 6 * h_square
        parity_allowed = c2 % 4 in (0, 3) if (plane or ramified_vertex) else c2 % 2 == 0
        if intersection >= spinor_marks and parity_allowed:
            result.append((n, intersection, c2))
    return result


assert candidates(1, 7, plane=True) == [(6, 6, 0), (7, 8, -1)]
assert candidates(2, 13) == [(12, 12, 0)]
assert candidates(2, 13, ramified_vertex=True) == [(12, 26, -1)]
xi, c1, c2, h = symbols("xi c1 c2 h")
polynomial = expand((2 * xi + h - c1)**3).as_poly(xi)
pushforward = polynomial.nth(3) * (c1**2 - c2) + polynomial.nth(2) * c1 + polynomial.nth(1)
assert expand(pushforward - 6 * h**2 - 2 * (c1**2 - 4 * c2)) == 0
assert 26 - (2 * 11 + 2) == 2
print("PASS: plane 6,7; smooth quadric/cone 12; ramified cover 24 distinct marks in length 26")
