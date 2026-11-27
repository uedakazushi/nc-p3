"""Slope of transitive cuts and distinct paths in the two-block tangent calculation."""
from itertools import combinations_with_replacement, product


def main():
    for dimension, length in product(range(1,6), (4,8,12)):
        for sub in combinations_with_replacement(range(dimension+1),length):
            if max(sub) == 0 or min(sub) == dimension:
                continue
            if any(0 < left == right < dimension for left, right in zip(sub, sub[1:])):
                continue
            assert sum(sub[:-1]) - (length-1)*sub[-1] < 0
    paths = {''.join(word) for length in range(10)
             for word in product(('x', 'yz'), repeat=length)}
    assert len(paths) == 2**10 - 1
    print('Stable cut inequalities and distinct tangent paths through nine concatenations verified.')


if __name__ == '__main__':
    main()
