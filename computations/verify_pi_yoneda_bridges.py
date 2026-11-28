"""Enumerate the vertex Euler and Yoneda dimension inequalities."""
from itertools import product
from math import comb


def feasible(d):
    local_euler = [2*d[i]-4*d[i-1]-4*d[(i+1)%4]+6*d[(i+2)%4] for i in range(4)]
    square = sum(x*y for x,y in zip(d, local_euler))
    minimum = 0
    for i in range(4):
        max_rows = 3*(d[i-1]-1)-d[i]
        max_cols = 3*(d[(i+1)%4]-1)-d[i]
        products = [a*b for a in range(max_rows+1) for b in range(max_cols+1)
                    if a+b+local_euler[i] >= 0]
        if not products:
            return False
        minimum += min(products)
    return minimum <= square


def main():
    surviving = set()
    for d in product((2,3), repeat=4):
        if len(set(d)) > 1 and feasible(d):
            surviving.add(min(d[j:]+d[:j] for j in range(4)))
    assert surviving == {(2,2,2,3)}
    for a,b in product(range(1,4), repeat=2):
        assert max(2, abs(a-b)+1) < min(comb(l+3,3) for l in range(1,4))
    print("All fourteen mixed profiles and all two-end obstruction matrix sizes checked.")


if __name__ == "__main__":
    main()
