"""Verify the binary center product and interval count used for rank nine."""
from itertools import product

from sympy import Matrix, symbols


def main():
    # With bases e,f and u,v the last slices are eu, ev+fu, fv.
    binary = Matrix([[1, 0, 0], [0, 1, 0], [0, 1, 0], [0, 0, 1]])
    assert binary.rank() == 3
    quotient = binary.T.nullspace()[0].reshape(2, 2)
    assert quotient.det() != 0
    s, t = symbols("s t")
    pencil = Matrix([[s, t, 0], [0, s, t]])
    minors = [pencil[:, [i, j]].det() for i, j in ((0, 1), (0, 2), (1, 2))]
    assert minors == [s**2, s*t, t**2]
    for length in range(2, 10):
        for ends in product((0, 1), repeat=length-1):
            assigned = [i+side for i, side in enumerate(ends)]
            if len(set(assigned)) == len(assigned):
                assert len(set(range(length))-set(assigned)) == 1
    print("PASS: binary center slices, perfect binary product, distinct plane assignments")


if __name__ == "__main__":
    main()
