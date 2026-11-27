"""Possible kernel splittings of a ten-column evaluation on the binary cubic."""
from itertools import combinations_with_replacement


def splittings(rank):
    count = 10-rank
    for degree in range(count, 3*rank+1):
        for parts in combinations_with_replacement(range(1, degree-count+2), count):
            if sum(parts) == degree:
                yield degree, parts


def main():
    assert not list(splittings(1)) and not list(splittings(2))
    for degree, parts in splittings(3):
        assert max(parts) <= 3
        for torsion in range(10-degree):
            free_degree = 12-degree-torsion
            assert free_degree >= 3 and free_degree-2+torsion <= 3
    possibilities = []
    for degree, parts in splittings(4):
        torsion = 12-degree
        fifth = torsion+sum(max(d-3, 0) for d in parts)
        if fifth >= 4:
            possibilities.append(torsion)
            if torsion <= 3:
                assert parts == (1, 1, 1, 1, 1, 7-torsion) and fifth == 4
            else:
                assert fifth == torsion and max(parts) <= 3
    assert set(possibilities) == set(range(7))
    print("All evaluation ranks, splitting degrees and fifth-moment dimensions verified.")


if __name__ == "__main__":
    main()
