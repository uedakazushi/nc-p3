"""Balanced rectangular sections and the coordinate-line diagrams of a pure frame."""
from collections import defaultdict
from itertools import combinations
from math import comb


def main():
    for a in range(2, 5):
        for b in range(2, 5):
            cuts = a+b-3
            assert (comb(cuts, a-2) == comb(cuts, a-1)) == (a == b)
    edges = list(combinations(range(4), 2))
    possibilities = defaultdict(set)
    for count in range(1, len(edges)+1):
        for subset in combinations(edges, count):
            degrees = tuple(sorted(sum(v in edge for edge in subset) for v in range(4)))
            if max(degrees) <= 2 and max(count-d for d in degrees) <= 2:
                possibilities[count].add(degrees)
    assert dict(possibilities) == {
        1: {(0, 0, 1, 1)}, 2: {(0, 1, 1, 2), (1, 1, 1, 1)},
        3: {(1, 1, 2, 2)}, 4: {(2, 2, 2, 2)},
    }
    print("Balanced section degrees and all five nonempty four-vertex diagrams verified.")


if __name__ == "__main__":
    main()
