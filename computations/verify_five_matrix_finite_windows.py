"""Rank and Chern bounds used to pass from five-matrix to four-matrix identities."""

import sympy as sp


def main():
    # Finite-edge restriction dimensions for supports of dimensions 0 through 4.
    supported = (0, 0, 1, 3, 6)
    gains = {}
    for rank in (2, 3):
        compressed = 6 - supported[4-rank]
        gains[rank] = 6 + compressed - supported[rank]
    assert gains == {2: 10, 3: 9}
    assert 16 - min(gains.values()) - 4 == 3

    h, c = sp.symbols("h c")
    inverse = sp.series(1/(1+h+c*h**2), h, 0, 4).removeO().expand()
    assert inverse.coeff(h, 3) == 2*c-1
    assert sp.expand((1+h)**4).coeff(h, 3) == 4

    bad_sizes = {size for size in (2, 3, 4) if 1+size+1 <= 5}
    assert bad_sizes == {2, 3}
    assert {size for size in bad_sizes if 2+size+1 <= 5} == {2}
    assert all(3+size+1 > 5 for size in bad_sizes)
    print("Four-dimensional extension ranks and five-matrix middle-block bounds verified.")


if __name__ == "__main__":
    main()
