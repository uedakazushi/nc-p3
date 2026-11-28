"""Conductor lengths on the cubic scroll and the four-base smoothing."""
from verify_cofactor_nineteen_clusters import intersection


def genera(m):
    square = 12*m-3*m*m
    canonical_degree = -m-2*(6-m)
    return 1+(square+canonical_degree)//2


assert [(genera(m), 6-m) for m in (1, 2, 3)] == [(0, 5), (2, 4), (1, 3)]
assert [m for m in (1, 2, 3) if 2*(3-genera(m)) <= 6-m] == [2]
for weights in [(0, 19, 293), (0, 23, 317)]:
    points = [(*(w-v for w in weights if w != v), -4*v) for v in weights]
    for multiplicity in (2, 1, 1, 1, 1):
        u, v, h = points.pop(0)
        points.extend([(u, v-u, h-multiplicity*u), (u-v, v, h-multiplicity*v)])
    assert intersection(points) == 16
print("PASS: conductor pair degrees and triple integral on the cubic-scroll smoothing")
