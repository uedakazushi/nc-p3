#!/usr/bin/env python3
"""Enumerate equality patterns of the three rank-two planes at each of four positions."""
import argparse
from collections import Counter
from itertools import combinations, permutations, product
import json
from pathlib import Path

import sympy as s


def pattern_counts():
    positions = range(4)
    neighbors = [[j for j in positions if j != i] for i in positions]
    # Canonical labels for all five set partitions of a three-element set.
    partitions = ((0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (0, 1, 2))
    counts = Counter()
    for rows in product(partitions, repeat=4):
        planes = [dict(zip(indices, row)) for indices, row in zip(neighbors, rows)]
        for i, j, k in permutations(positions, 3):
            if planes[k][i] != planes[k][j] and planes[i][j] != planes[i][k]:
                break
        else:
            shared = sum(len(set(row)) == 1 for row in rows)
            counts[shared] += 1
    return dict(sorted(counts.items()))


def binary_pencil_minors():
    # The two rows are the coefficients of x and y after contracting the source.
    a, b = s.symbols("a b")
    pencil = s.Matrix([[a, b, 0], [0, a, b]])
    minors = [pencil[:, [i, j]].det() for i, j in combinations(range(3), 2)]
    assert minors == [a**2, a*b, b**2]
    return [str(minor) for minor in minors]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    counts = pattern_counts()
    assert counts == {2: 12, 3: 16, 4: 1}
    result = {"patterns_tested": 5**4, "compatible_patterns": sum(counts.values()),
              "by_number_of_common_kernel_positions": counts,
              "projected_binary_pencil_minors": binary_pencil_minors()}
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
