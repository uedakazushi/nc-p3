"""Cyclic quartic eigenlines and the infinitesimal Segre scheme of the multiplication space."""
from itertools import combinations_with_replacement

import sympy as sp


def main():
    x, e = sp.symbols("x e")
    for count in range(1, 5):
        for parts in combinations_with_replacement(range(1, 5), count):
            if sum(parts) != 4:
                continue
            blocks = []
            for root, size in enumerate(parts):
                jordan = sp.eye(size)*root + sp.Matrix(size, size, lambda i, j: int(j == i+1))
                blocks.append(jordan)
                eigen = sp.Matrix([e**i for i in range(size)])
                assert all(sp.rem(z, e**size, e) == 0 for z in (jordan-(root+e)*sp.eye(size))*eigen)
                # The dual multiplication matrices are constant on antidiagonals.
                value = eigen*eigen.T
                for i in range(size):
                    for j in range(size):
                        assert sp.rem(value[i, j], e**size, e) == (e**(i+j) if i+j < size else 0)
            pencil = x*sp.eye(4)-sp.diag(*blocks)
            for root in range(count):
                assert pencil.subs(x, root).rank() == 3
    print("All length-four cyclic pencils, including the local quadruple, retain their eigenline and Segre jets.")


if __name__ == "__main__":
    main()
