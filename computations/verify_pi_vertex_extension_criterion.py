"""Euler and rank conditions for the four vertex Koszul complexes."""
from itertools import product
from math import comb


def euler_vector(dimensions):
    return tuple(sum((-1)**j * comb(4, j) * dimensions[(i+j) % 4]
                     for j in range(5)) for i in range(4))


def main():
    for dimensions in product(range(5), repeat=4):
        euler = euler_vector(dimensions)
        assert sum(euler) == 0
        assert (euler == (0,)*4) == (len(set(dimensions)) == 1)
    for dimension in range(1, 5):
        ranks = [dimension, 3*dimension, 3*dimension, dimension]
        terms = [comb(4,j)*dimension for j in range(5)]
        assert all(terms[j] == left+right for j,(left,right) in
                   enumerate(zip([0]+ranks, ranks+[0])))
    print('Vertex Euler sums, balanced dimension vectors, and exact-complex ranks checked.')


if __name__ == '__main__':
    main()
