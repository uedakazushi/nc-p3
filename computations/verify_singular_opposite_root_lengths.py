"""Check the small Hilbert-function and local-length bounds for singular opposite roots."""

from itertools import combinations, combinations_with_replacement


def shadow(words, variables):
    return {tuple(sorted(word + (i,))) for i in range(variables) for word in words}


def main():
    ternary = list(combinations_with_replacement(range(3), 2))
    assert min(len(shadow(basis, 3)) for basis in combinations(ternary, 4)) == 8
    # Spend n-1 units of total nonreduced length at each planar local factor.
    max_defect = [0]
    for excess in range(1, 8):
        max_defect.append(max(max_defect[excess-j]+j//2 for j in range(1, excess+1)))
    assert max_defect == [0, 0, 1, 1, 2, 2, 3, 3]
    # Curvilinear condition summed over Sing(h), of length at most three.
    assert min((20-d+1)//2 for d in range(4)) == 9
    print("Minimum ternary shadow 8; planar closure defects at excess up to seven bounded by 3.")


if __name__ == "__main__":
    main()
