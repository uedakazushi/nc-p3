"""Check the compressed Hilbert multiplicity and the full-cycle weight decomposition."""

from math import comb

import sympy as sp


def main():
    t = sp.symbols("t")
    numerator = sum(sum((-1)**j*comb(4, j)*comb(4*(i-j)+3, 3)
                        for j in range(i+1))*t**i for i in range(4))
    assert numerator.subs(t, 1) == 64
    for order in (3, 5, 7, 11):
        scalars = list(range(1, order+1))
        matrix = sp.zeros(order)
        for j, scalar in enumerate(scalars):
            matrix[(j+1) % order, j] = scalar
        assert matrix**order == sp.prod(scalars)*sp.eye(order)
    print("Compressed multiplicity 64 and scalar-arrow cyclic monodromy verified.")


if __name__ == "__main__":
    main()
