"""Check the linear restriction class, factor partitions and marked multiplicities."""
from math import comb

def partitions(total, maximum):
    if total == 0:
        yield ()
    for degree in range(1, min(maximum, total) + 1):
        for tail in partitions(total - degree, degree):
            yield (degree, *tail)

length = comb(4 + 3 - 1, 3)
assert length == 20
for marks in range(18, 21):
    assert length - marks < 5 - 2
    quartic_excess = max(16 + 3, 17 + 1 + 1, 3 + 3 + 1)
    bounds = {1: 7, 2: 15, 3: 24, 4: marks + quartic_excess}
    cases = list(partitions(6, 4))
    assert len(cases) == 9
    assert all(sum(bounds[d] for d in p) < 3*marks for p in cases)
assert -3 + 2 < 0
print("PASS: restriction length 20; all nine low-degree partitions violate the triple-mark bound")
