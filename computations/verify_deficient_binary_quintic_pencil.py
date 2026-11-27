"""Two-site cuts on all cyclic torsion divisors of length two and three."""
from itertools import combinations_with_replacement, product

import sympy as sp


def main():
    for length in (2, 3):
        for count in range(1, length+1):
            for parts in combinations_with_replacement(range(1, length+1), count):
                if sum(parts) != length:
                    continue
                for vanishing in product(*(range(n+1) for n in parts)):
                    if sum(vanishing) != 2:
                        continue
                    blocks = [sp.Matrix(n, n, lambda i, j: int(i == j+order))
                              for n, order in zip(parts, vanishing)]
                    rank = sp.diag(*blocks).rank()
                    assert rank == length-2
                    assert rank+4-length == 2
    print("Both endpoint rank bounds verified on every divisor of length two in length two or three.")


if __name__ == "__main__":
    main()
