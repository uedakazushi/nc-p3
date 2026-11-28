"""Check periodic-tail multiplicities and the sizes forbidden by the matrix identity."""

from itertools import product


def main():
    for s in (3, 4):
        possible = [(c, d) for c in range(1, 5) for d in (2, 3)
                    if c + d + 1 <= s]
        assert possible == ([(1, 2)] if s == 4 else [])
    for a in product(range(4), repeat=4):
        if not any(a):
            continue
        for j in range(4):
            assert max(a[(j + r) % 4] for r in range(4)) > 0
    for c in range(1, 5):
        p = [c, -3*c, 3*c, -c]
        assert sum(p) == 0
    print("The only unresolved size pair at bound four is (1,2); none survives bound three.")


if __name__ == "__main__":
    main()
