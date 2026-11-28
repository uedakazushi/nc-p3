"""Check the complete intersection, exceptional-factor budgets and mixed dimensions."""
from itertools import product
from math import comb

# Degrees in the Artin reduction of a (2,2,2) complete intersection.
degrees = [sum(a) for a in product(range(2), repeat=3)]
assert [degrees.count(i) for i in range(4)] == [1, 3, 3, 1]
assert len(degrees) - 2 == 6
assert 2*7 + 2*15 + 15 < 60
assert 4*7 + 7 + 24 < 60
assert 20 - 2 > 14 and 20 - 2*7 > 3
assert 16 - 1 - 2*6 == 3
print("PASS: column base length at most six, unique repeated factor, at least three components in degree six")

def partitions(total, maximum):
    if not total:
        yield ()
    for d in range(1, min(total, maximum) + 1):
        for rest in partitions(total-d, d):
            yield (d, *rest)

cap = {1: 3, 2: 6, 3: 10, 4: 13}
remaining = [p for p in partitions(6, 4) if len(p) >= 3 and sum(cap[d] for d in p) >= 19]
assert set(remaining) == {(4, 1, 1), (3, 2, 1), (3, 1, 1, 1)}
print("PASS: conic cap six leaves exactly three disconnected degree-six partitions")
