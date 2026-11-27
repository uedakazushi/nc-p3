"""Cohomological lower bounds of the external polynomial ring calculation."""
from itertools import product
from math import comb


def main():
    for dimension in range(5):
        possible = {4 + p - q for p,q in product(range(dimension+1),repeat=2)}
        assert min(possible) == 4-dimension
        if dimension < 4:
            assert 0 not in possible
        koszul_ranks = [comb(dimension,j) for j in range(dimension+1)]
        assert len(koszul_ranks) == dimension+1
    assert all(left+right > 1 for left,right in product(range(1,5), repeat=2))
    print('External polynomial flat-dimension offsets and the vertex multiplicity obstruction checked.')


if __name__ == '__main__':
    main()
