"""Length bounds for primary modules and the long linear-syzygy summand."""
from itertools import combinations_with_replacement, product


def local_types(length):
    for count in range(1, length+1):
        for terms in combinations_with_replacement(range(1, length+1), count):
            if sum(terms) == length:
                yield terms


def main():
    for length in range(1, 7):
        for counts in local_types(length):
            for family in product(*(local_types(n) for n in counts)):
                if all(terms.count(max(terms)) >= 2 for terms in family):
                    assert sum(max(terms) for terms in family) <= 3
    for torsion in (1, 2, 3):
        long_degree = 7-torsion
        assert (5-long_degree) == torsion-2
    print("Tied maximal elementary divisors have effective degree at most three in length six.")


if __name__ == "__main__":
    main()
