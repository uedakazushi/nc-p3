"""Curve span bounds and a finite complementary pencil at a conic of pure relations."""
from itertools import product

import sympy as sp


def main():
    pairs = [(a, b) for a, b in product(range(2, 5), repeat=2)
             if a+b-1 <= 6 and (a == b == 4 or a+b-1 <= 4)
             and (min(a, b) != 2 or a+b-1 <= 3)]
    assert pairs == [(2, 2)]
    u, v = sp.symbols("u v")
    pencil = sp.diag(u, v, u-v).col_join(sp.zeros(1, 3))
    assert sp.expand(pencil[:3, :].det()-u*v*(u-v)) == 0
    assert all(pencil.subs({u: a, v: b}).rank() == 2 for a, b in [(1, 0), (0, 1), (1, 1)])
    print("Only the (2,2) moving support survives; the complementary pencil has cubic maximal minors.")


if __name__ == "__main__":
    main()
