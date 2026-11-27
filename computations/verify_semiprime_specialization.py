"""Sizes in the regular, generic exceptional, and special simple quotient calculations."""
from itertools import product


def main():
    matrix_sizes = [(m, q) for m, q in product(range(1, 5), repeat=2) if m*q == 4]
    assert matrix_sizes == [(1,4), (2,2), (4,1)]
    singular_profiles = {(1,3), (2,2)}
    assert all(sum(parts) == 4 for parts in singular_profiles)
    assert all(max(parts) < 4 for parts in singular_profiles)
    assert 4*48 - 4 == 188
    print("Degree-four simple quotients, smaller exceptional simple sizes, and the common parameter socle.")


if __name__ == "__main__":
    main()
