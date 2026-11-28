"""Local fixed-line thickenings used to exclude rank-two quadratic images."""

from itertools import product


def divides(a, b):
    return all(i <= j for i, j in zip(a, b))


def square(generators):
    return [tuple(i+j for i, j in zip(a, b)) for a, b in product(generators, repeat=2)]


def main():
    bounds = {n: (5-n)*length for n, length in [(1, 1), (2, 3), (3, 7)]}
    partitions = [(3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]
    assert [sum(bounds[i] for i in part) for part in partitions] == [18, 18, 17, 16]

    normal_square = square([(0, 1), (2, 0)])
    standard = [m for m in product(range(4), range(2))
                if not any(divides(n, m) for n in normal_square)]
    assert len(standard) == 6
    # In the (3,1) chart every exceptional monomial has a-degree <2, b-degree <6
    # and c-degree <2: a^2*b^2, a*b*c, b^6, c^2 belong to the ideal square.
    generators = [(0, 0, 1), (0, 3, 0), (1, 1, 0)]
    special = [m for m in product(range(2), range(6), range(2)) if m[1]+m[2] >= 2
               and not any(divides(n, m) for n in square(generators))]
    assert len(special) == 8
    assert 3*3+len(special) == 17
    print("Multiplicity bounds verified for the semisimple decompositions and nilpotent Jordan types.")


if __name__ == "__main__":
    main()
