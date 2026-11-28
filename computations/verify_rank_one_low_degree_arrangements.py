"""Check ruling degrees, common sections and the low-curve partition budgets."""
from itertools import product

assert 3*3 > 8 and 6-2*3 == 0 and 4*14 < 60
for d, components, shared in [(5, 1, 2), (4, 1, 3), (4, 2, 2)]:
    assert 8 - (d+components) == shared
assert 3+2+2 < 16-8
for lines_a, lines_b, conics in product(range(6), repeat=3):
    if lines_a+lines_b+2*conics != 5:
        continue
    a, b = sorted((lines_a+conics, lines_b+conics))
    if not a:
        continue
    if a == 1:
        assert max(lines_a, lines_b) >= 3
        assert max(lines_a, lines_b) >= 4 or conics == 1
    else:
        assert (a,b) == (2,3) and max(lines_a, lines_b) > 0
print("PASS: ruling intersections, common sections, four-slope bound and all quadric bidegrees")
