"""Check sharp polynomial path growth and a mixed-endpoint exponential obstruction."""
from math import comb

from sympy import Matrix, ones


def main():
    adjacency = Matrix(4, 4, lambda a, b: comb(b-a+3, 3) if a <= b else 0)
    for blocks in range(21):
        paths = (ones(1, 4)*adjacency**blocks*ones(4, 1))[0]
        assert paths == comb(4*blocks+4, 3)
    obstruction = Matrix([[1, 1], [1, 0]])
    # Returning in two steps already has two independent choices.
    assert (obstruction**2)[0, 0] == 2
    assert (obstruction**48)[0, 0] >= 2**24 > comb(4*48+4, 3)
    print("PASS: four chained cycles attain cubic growth; mixed-endpoint returns exceed it")


if __name__ == "__main__":
    main()
