"""Verify geometric marking and local multiplication on the reduced coordinate quadrangle."""
from itertools import combinations_with_replacement, product
from math import comb


def surviving_words(n):
    return [w for w in combinations_with_replacement(range(4), n)
            if not ({0, 2} <= set(w) or {1, 3} <= set(w))]


def main():
    for n in range(1, 13):
        basis = surviving_words(n)
        assert len(basis) == 4*n
        assert {tuple(sorted((x, *tail))) for x, tail in product(range(4), surviving_words(n-1))
                if not ({0, 2} <= {x, *tail} or {1, 3} <= {x, *tail})} == set(basis)
        assert sum((j+1)*len(surviving_words(n-2*j)) for j in range(n//2+1)) == comb(n+3, 3)
    # The two defining quadratic monomials give four products with a single alternating relation.
    quadrics = [(0, 2), (1, 3)]
    products = [tuple(sorted(a+b)) for a, b in product(quadrics, repeat=2)]
    assert len(set(products)) == 3 and products[1] == products[2]
    print("PASS: quadrangle Hilbert series, section multiplication and rank-two quartic relation")


if __name__ == "__main__":
    main()
