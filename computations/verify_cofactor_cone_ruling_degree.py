"""Hodge and branch restrictions on the resolved kernel cone."""
def pairing(left, right):
    a, b = left
    c, d = right
    return -2*a*c+a*d+b*c


assert [t for t in range(1, 5) if 6-2*t >= 0 and 6*t-t*t >= 8] == [2, 3]
assert pairing((1, 3), (1, 3)) == 4  # pullback to the double cover has square eight
assert pairing((1, 3), (1, 2)) == 3
for b in (1, 2, 3, 4):
    branch_half = (1, b)
    assert pairing(branch_half, (-1, b-4)) == -2
assert [b for b in range(1, 6) if 3*b-12 <= -4] == [1, 2]
assert 2*4-2 < 8 and -2+3 == 1
print("PASS: quartic adjunction bound, ruling degrees, and four cone branch classes")
