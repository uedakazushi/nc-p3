"""Check the two-component adjunction degrees and the remaining degree partitions."""
from itertools import combinations_with_replacement, product


def main():
    for genus in (0, 1):
        length = 5 + 1 - genus
        canonical = (-2 + length, 2*genus - 2 + length)
        assert (canonical[0] - 3, canonical[1] - 4) == (1-genus, genus)
        assert min(length - e for e in range(1, 5)) > 0
    for labels in product(range(3), repeat=4):
        i, j, k, l = labels
        if sorted((i, k, k)) == sorted((j, j, l)):
            assert k == j and i == l
    for total in (5, 6, 7):
        decompositions = {v for n in range(1, 4)
                          for v in combinations_with_replacement(range(3, total+1), n)
                          if sum(v) == total}
        assert decompositions == {5: {(5,)}, 6: {(6,), (3, 3)}, 7: {(7,), (3, 4)}}[total]
    print("PASS: two-component residual degrees, divisor recurrence, and degree partitions")


if __name__ == "__main__":
    main()
