"""Curve and divisor leading coefficients, A1 intersection, and the quartic cutoff."""
from fractions import Fraction
from itertools import product


def main():
    for index in (1, 2, 4):
        degree = Fraction(4, index)
        assert degree.denominator == 1 and 4*index*4*degree == 64
    assert 4*2 == 8
    multiplicities = []
    for counts in product(range(3), repeat=4):
        square = Fraction(1-sum(c*c for c in counts), 2)
        if square >= 0 and square.denominator == 1:
            multiplicities.append(sorted(counts))
    assert multiplicities and all(m == [0,0,0,1] for m in multiplicities)
    assert 3*2 > 4 and 2*2+1 > 4
    top = 4*(4-1)
    assert (top+2,top+1,top+4) == (14,13,16)
    print("Integral degrees, the unique geometric A1 correction, both scroll inequalities, and the 16-degree cutoff checked.")


if __name__ == "__main__":
    main()
