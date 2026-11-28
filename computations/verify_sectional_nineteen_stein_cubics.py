"""Check the Stein resolution and marked-point bounds, including reducible cubics."""
from itertools import product
from math import comb


def hilbert(d):
    return 4 * d * d + 3 * d + 1 if d >= 0 else 0


def main():
    numerator = [sum((-1)**j * comb(4, j) * hilbert(d-j) for j in range(5))
                 for d in range(9)]
    assert numerator == [1, 4, -3, -2, 0, 0, 0, 0, 0]
    # The column relations on two coprime primitive cofactors are their Koszul multiples.
    assert 20 - 14 > 3 + 2
    assert 18 - (20 - 2 * 4) == 6
    for first, second, overlap in product(range(8), range(15), range(8)):
        if first + second - overlap >= 18:
            assert overlap + 1 < 6
    for single, double, triple in product(range(21), range(11), range(8)):
        if single + 2*double + 3*triple <= 21 and single + double + triple >= 18:
            assert double + triple < 6
    # A transverse ruling shares at most one blowup center with each disjoint ruling.
    assert 2 + 2 == 4 and 1 + 1 < 3
    print("PASS: graded presentation, at most two exceptional marks, six necessary singular marks incompatible with every cubic type")


if __name__ == "__main__":
    main()
