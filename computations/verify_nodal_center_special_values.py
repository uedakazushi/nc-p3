"""Local rank partitions and Euler/Chern numbers for both central polarizations."""
from itertools import combinations_with_replacement, product
from math import gcd

import sympy as sp


def main():
    pairs = list(product(range(1, 5), repeat=2))
    for size in range(1, 5):
        for rows in combinations_with_replacement(pairs, size):
            if sum(d*m for d, m in rows) == 4 and gcd(*(m for d, m in rows)) == 1:
                assert min(m for d, m in rows) == 1
    h = sp.symbols("h")
    c = sp.series((1+h)**5 / ((1+2*h)*(1+4*h)), h, 0, 3).removeO()
    smooth_euler = 2*4 - 8*c.coeff(h, 2)
    nodes = (4-smooth_euler)//2
    blowup_euler = smooth_euler + 4*nodes
    assert (smooth_euler, nodes, blowup_euler, (blowup_euler-2)//2) == (-56, 30, 64, 31)
    print("Primitive rank bound, both Euler counts, thirty nodes and blowup Betti number 31.")


if __name__ == "__main__":
    main()
