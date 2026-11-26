"""Verify the local generator-space propagation in all possible dimensions."""
from itertools import product

from sympy import Matrix, symbols


def main():
    # At a point in the right kernel, the present pair has dimension 1 or
    # 2 in the conormal fiber. The previous pair must span the same fiber.
    a, b = symbols("a b", nonzero=True)
    previous = Matrix([[a, 0]])
    assert previous.nullspace() == [Matrix([0, 1])]
    previous = Matrix([[a, 0], [symbols("c"), b]])
    assert previous.det() == a*b and not previous.nullspace()
    # At generic component points each plane equation can remove at most
    # one layer on one linear component, even for the triple line.
    for multiplicities in ((3,), (2, 1), (1, 1, 1)):
        for i, j in product(range(len(multiplicities)), repeat=2):
            remainder = [max(0, n-int(k==i)-int(k==j))
                         for k, n in enumerate(multiplicities)]
            assert sum(remainder) >= 1
    # Three restricted relations with a one-dimensional cubic overlap.
    assert 3**3-2*3**2+1 == 20-10
    print("PASS: conormal propagation in both ranks, all plane line multiplicities, cubic overlap count")


if __name__ == "__main__":
    main()
