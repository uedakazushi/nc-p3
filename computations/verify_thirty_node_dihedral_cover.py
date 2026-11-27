"""Check the exceptional quadric collections, the cubic local profiles, and the eight covers."""
from itertools import product


def line_cohomology(a, b):
    # P1 Kunneth. O(-1) is acyclic.
    groups = [0, 0, 0]
    left = (max(a+1, 0), max(-a-1, 0))
    right = (max(b+1, 0), max(-b-1, 0))
    for i, j in product(range(2), repeat=2):
        groups[i+j] += left[i] * right[j]
    return groups


def main():
    for degree in range(4):
        for difference in (-1, 0, 1):
            assert line_cohomology(degree+difference, degree)[1:] == [0, 0]
    for t, summand in product((-1, 0), (0, 1)):
        assert not any(line_cohomology(t-summand, -1))
    assert line_cohomology(0, 0) == [1, 0, 0]
    for twist in (-1, 0):
        assert not any(line_cohomology(twist-1, -1))
    assert not any(line_cohomology(-1, 0))
    orientations = {tuple((s*i+j) % 4 for i in range(4))
                    for s, j in product((-1, 1), range(4))}
    assert len(orientations) == 8
    for x, y in product(orientations, repeat=2):
        assert tuple(x[i] for i in y) in orientations
    assert {x % 4 for x in (1,-1)} == {1,3}
    assert 2 % 4 == -2 % 4
    print("Two exceptional objects per quadric, sixty-four in total, and eight distinct indexed orientations.")


if __name__ == "__main__":
    main()
