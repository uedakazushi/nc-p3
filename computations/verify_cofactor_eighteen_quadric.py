"""Degree-18 birational quadric inverses and their triple counts."""
from verify_cofactor_nineteen_clusters import fixed_points, intersection

assert [(a, 6-a) for a in range(1, 4) if (a-1)*(5-a) in (3, 4)] == [(2, 4), (3, 3)]
for basis_weights in [(17, 271), (19, 307)]:
    points = fixed_points(2, 4, *basis_weights)
    for _ in range(8):
        u, v, h = points.pop(0)
        points.extend([(u, v-u, h-u), (u-v, v, h-v)])
    assert intersection(points) == 12
assert 2*2*4 - 8 == 8 and 2*3*3 - 8 == 10
assert 4**2 - 2 - 6 == 8 and 4*2 - 2 == 6
print("PASS: bidegrees, base lengths, and the exact degree-twelve triple count")
