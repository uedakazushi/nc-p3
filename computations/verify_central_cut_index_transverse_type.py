"""Integral Hilbert multiplicities and transverse representation dimensions."""
from itertools import combinations_with_replacement, product

import sympy as sp


def main():
    m, i = sp.symbols("m i")
    polynomial = (4*m+i+1)*(4*m+i+2)*(4*m+i+3)/6
    assert sp.expand(polynomial).coeff(m, 3)*6 == 64
    possibilities = [(degree, rank) for degree, rank in product(range(1, 65), range(1, 9))
                     if degree*rank**2 == 64]
    assert possibilities == [(1,8), (4,4), (16,2), (64,1)]
    for degree, rank in possibilities:
        assert [r for r in range(1, rank+1) if degree*rank*r % 64 == 0] == [rank]
    rows = list(product(range(2, 5), range(1, 5)))
    solutions = [v for n in range(2, 5) for v in combinations_with_replacement(rows, n)
                 if sum(d*r for d, r in v) == 4]
    assert solutions == [((2,1),(2,1))]
    print("Four normalization degrees, full division indices, and unique nontrivial transverse representation profile.")


if __name__ == "__main__":
    main()
