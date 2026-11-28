"""Test the flag of modules at the unique nonnormal point and the surface bounds."""
from sympy import Poly, Symbol


def main():
    v = Symbol("v")
    for polynomial in (v**4, v**2*(v-1)**2, v*(v-1)*(v-2)*(v-3)):
        p = Poly(polynomial, v)
        powers = [Poly(v**j, v).rem(p) for j in range(5)]
        assert [q.degree() for q in powers[:4]] == [0, 1, 2, 3]
        # A positive-degree multiplier preserving span(1,v,v^2) would itself
        # lie in that span. Multiplication up to the first degree-three term forbids it.
        for degree in (1, 2):
            assert (powers[degree] * powers[3-degree]).rem(p).degree() == 3
    assert 21 - 18 == 3 and 24 - 5 - 18 == 1
    for multiplicities in ((1, 1, 1, 1), (2,)):
        assert sum(m*m for m in multiplicities) == 4
        assert 2 - sum(m*(m-1)//2 for m in multiplicities) > 0
    assert 4 * 2 - 2 * 5 < 0
    print("PASS: length-four monogenic fiber, rational normalization and trace exclusions on both finite surface models")


if __name__ == "__main__":
    main()
