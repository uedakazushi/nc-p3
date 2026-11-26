"""Enumerate the remaining corridor lengths with the two crossing capacities."""
from collections import Counter
from itertools import product


def fits(lengths):
    ends = {(origin + n + 1) % 4 for origin, n in zip((0, 2), lengths)}
    if ends != {0, 2}:
        return False
    moving_edges = Counter((origin + j) % 4 for origin, n in zip((0, 2), lengths)
                           for j in range(1, n))
    return all(moving_edges[i] <= 2 for i in (1, 3))


allowed = {tuple(sorted(pair)) for pair in product(range(3, 41), repeat=2) if fits(pair)}
assert allowed == {(3, 3), (3, 7), (5, 5)}
print("PASS: the only unordered lengths are (3,3), (3,7), (5,5)")
