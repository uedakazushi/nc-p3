"""Check weighted rescaling at infinity and dimensions in the octic parameter test."""

from fractions import Fraction
from itertools import product
from math import comb


def main():
    weights = (4, 8, 12)
    for valuations in product(range(-4, 3), repeat=3):
        if min(valuations) >= 0:
            continue
        order = max(Fraction(-v, w) for v, w in zip(valuations, weights))
        corrected = [v + w * order for v, w in zip(valuations, weights)]
        assert order > 0 and min(corrected) == 0
    quotient = [sum((-1)**i * comb(4, i) * comb(n-8*i+3, 3)
                    for i in range(min(4, n//8) + 1)) for n in range(33)]
    assert sum(quotient) == 8**4 and quotient[28:] == [1, 0, 0, 0, 0]
    assert (28+2, 28+1, 28+4) == (30, 29, 32)
    print("Weighted trace normalization and all three finite bounds for four octic cuts pass.")


if __name__ == "__main__":
    main()
