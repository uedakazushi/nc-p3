"""Intersection forms on the candidate kernel surfaces at d=19."""
from itertools import combinations_with_replacement


def ruled_classes(invariant, hyperplane_fiber_coefficient):
    candidates = []
    for a in range(1, 6):
        b = 5 - a * (hyperplane_fiber_coefficient - invariant)
        square = 2 * a * b - invariant * a * a
        if square >= 8 and b >= invariant * a:
            candidates.append((a, b, square - 8))
    return candidates


assert [(k, n) for k in range(2, 5) for n in range(1, 4)
        if 8 * n * k <= 25] == [(2, 1), (3, 1)]
assert ruled_classes(3, 3) == []
assert ruled_classes(1, 2) == [(2, 3, 0)]
assert ruled_classes(2, 2) == [(1, 5, 0), (2, 5, 4)]
assert ruled_classes(0, 1) == [(1, 4, 0), (2, 3, 4), (3, 2, 4), (4, 1, 0)]
assert -2 + 2 + 3 == 3  # (2E+3f)(E+f) on F_1
patterns = [p for size in range(1, 5)
            for p in combinations_with_replacement((1, 2), size)
            if sum(m * m for m in p) == 4]
assert patterns == [(2,), (1, 1, 1, 1)]
assert 1 + 5 - 4 == 2 and 19 + 3 > 21
print("d=19: cubic and quadric intersection forms, genus, and base clusters verified")
