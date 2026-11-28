"""Check the integer arguments in the stabilizer and cutoff proofs on small inputs."""

from itertools import combinations, combinations_with_replacement
from math import gcd

import sympy as sp


def check_common_degrees():
    for bound in range(1, 12):
        for size in range(1, bound + 1):
            for generators in combinations(range(1, bound + 1), size):
                divisor = gcd(*generators)
                reduced = [degree // divisor for degree in generators]
                first = reduced[0]
                threshold = (bound * bound + divisor - 1) // divisor
                reached = {0}
                for value in range(1, threshold + first):
                    if any(value - degree in reached for degree in reduced):
                        reached.add(value)
                assert all(value in reached for value in range(threshold, threshold + first))


def check_valuations():
    for length in range(1, 9):
        for dims in combinations_with_replacement(range(4), length):
            fixed = [*reversed(dims), 0]
            valuations = [next(i for i, value in enumerate(fixed) if value <= 4 - k)
                          for k in range(1, 5)]
            assert sum(valuations) == sum(dims)
    n = sp.symbols("n")
    hilbert = sp.expand_func(sp.binomial(4 * n + 3, 3)).as_poly(n)
    assert 6 * hilbert.LC() == 64


if __name__ == "__main__":
    check_common_degrees()
    check_valuations()
    print("Verified the common-degree bound, fixed-stratum sum, and Hilbert leading term.")
