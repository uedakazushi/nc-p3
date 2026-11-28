"""Length-three tautological intersection on four blowups in degree 19."""
from fractions import Fraction
from itertools import product

PARTITIONS = [(), (1,), (1, 1), (2,), (1, 1, 1), (2, 1), (3,)]


def local_terms(point, size):
    u, v, h = point
    for partition in PARTITIONS:
        if sum(partition) != size:
            continue
        weights = []
        denominator = 1
        for i, row in enumerate(partition):
            for j in range(row):
                arm = row - j - 1
                leg = sum(length > j for length in partition[i + 1:])
                weights.append(h - i * u - j * v)
                denominator *= ((arm + 1) * v - leg * u) * ((leg + 1) * u - arm * v)
        yield weights, denominator


def intersection(fixed):
    def integrate(index, remaining, weights, denominator):
        if index == len(fixed):
            if remaining:
                return 0
            a, b, c = weights
            c1, c2, c3 = a + b + c, a * b + a * c + b * c, a * b * c
            return Fraction(c2**3 - 2 * c1 * c2 * c3 + c3**2, denominator)
        return sum(integrate(index + 1, remaining - k, weights + w, denominator * d)
                   for k in range(remaining + 1)
                   for w, d in local_terms(fixed[index], k))
    return integrate(0, 3, [], 1)


def fixed_points(a, b, u, v):
    return [(u * (1 - 2 * i), v * (1 - 2 * j), -a * i * u - b * j * v)
            for i, j in product(range(2), repeat=2)]


def blow_up_corners(points):
    charts = []
    for u, v, h in points:
        charts.extend([(u, v - u, h - u), (u - v, v, h - v)])
    return charts


for parameters in [(17, 71), (13, 103)]:
    assert intersection(fixed_points(2, 2, *parameters)) == 20
    assert intersection(blow_up_corners(fixed_points(2, 3, *parameters))) == 16

print("PASS: exact length-three Hilbert-scheme intersections 20 and 16")
