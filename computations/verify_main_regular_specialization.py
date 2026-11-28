"""Enumerate reduced-trace fiber multiplicities and compute the quartic bounds."""

from itertools import product
from math import comb


def main():
    mixed = []
    for m in product(range(5), range(3), range(2), range(2)):
        d = sum(j * x for j, x in enumerate(m, 1))
        if not 0 < d < 4:
            continue
        assert m[3] == 0 and sum(m[1:]) <= 1
        mixed.append((m, 4 - d))
    assert len(mixed) == 6
    dim = [sum((-1)**j * comb(4, j) * comb(n-4*j+3, 3)
               for j in range(min(n//4, 4)+1)) for n in range(17)]
    assert dim == [1,4,10,20,31,40,44,40,31,20,10,4,1,0,0,0,0]
    last = max(n for n, h in enumerate(dim) if h)
    assert (last + 2, last + 1, last + 4) == (14, 13, 16)
    print("Six mixed semisimple multiplicity profiles; each has at most one nonthin factor.")
    print("Quartic parameter quotient and path, quotient, Koszul bounds verified.")


if __name__ == "__main__":
    main()
