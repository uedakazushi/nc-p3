"""Check the trace obstruction and rank-weighted sectional degree bounds."""
from itertools import product


def main():
    assert 2*4-6 == 2 and 6-4+1 == 3
    assert [degree for degree in range(2,8) if 8*degree <= 4**2] == [2]
    assert 6*8-2*5 == 38
    for rank_two,rank_one in product(range(22),range(8)):
        if 2*rank_two+6*rank_one <= 38:
            assert rank_two+3*rank_one <= 19
    for genus in range(5):
        for surplus in range(1,5):
            degree = genus+3+surplus
            assert 2*degree-(2*genus+6) == 2*surplus
            assert (2*genus+6-degree)+1-genus == 4-surplus
    print("PASS: trace divisor length and evaluation bounds, quartic Hodge equality, weighted degree at most nineteen")


if __name__ == "__main__":
    main()
