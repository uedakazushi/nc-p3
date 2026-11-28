"""Check the numerical alternatives in the sectional genus-three graph argument."""
from fractions import Fraction


def intersection(left, right):
    return left[0]*right[1] + left[1]*right[0]


def main():
    assert [(a, 6-a) for a in range(1, 4) if a*(6-a) >= 8] == [(2, 4), (3, 3)]
    h, l, canonical = (2, 2), (2, 1), (-2, -2)
    assert (intersection(h, h), intersection(h, l), intersection(l, l)) == (8, 6, 4)
    assert tuple(2*x-y-2*z for x, y, z in zip(l, canonical, h)) == (2, 0)
    assert 1 + (intersection(h, h) + intersection(canonical, h))//2 == 1
    assert 1 + (8-6)//2 == 2 < 3
    assert 8 - 2*Fraction(3, 2)*6 + Fraction(9, 4)*4 == -1
    assert 1 < 2 < 2**2
    assert max(2*g - 2 for g in (3, 4)) < 8
    print("PASS: both double-cover types, sectional genera, Hodge bound and quartic adjunction obstruction")


if __name__ == "__main__":
    main()
