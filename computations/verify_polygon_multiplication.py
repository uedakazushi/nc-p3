"""Rank lower bounds and complementary-chain certificates on common polygons."""
from collections import Counter
from itertools import combinations_with_replacement

from verify_binary_product_obstruction import distributions


def orbit(d, e):
    representatives = []
    for a, b in ((d,e), (e,d), (d[::-1],e[::-1]), (e[::-1],d[::-1])):
        representatives.extend((a[i:]+a[:i], b[i:]+b[:i]) for i in range(len(a)))
    return min(representatives)


def monomial_labels(degrees):
    """With gluing scalar one, match the coefficients on the two sides of every node."""
    n = len(degrees)
    nodes = list(range(n))

    def root(i):
        while nodes[i] != i:
            i = nodes[i]
        return i

    for i, degree in enumerate(degrees):
        if not degree:
            nodes[root(i)] = root((i+1) % n)
    classes = sorted({root(i) for i in range(n)})
    next_index = len(classes)
    labels = []
    for i, degree in enumerate(degrees):
        row = [classes.index(root(i))]
        if degree:
            row += list(range(next_index, next_index+degree-1))
            row.append(classes.index(root((i+1) % n)))
            next_index += degree-1
        labels.append(row)
    assert next_index == 4
    return labels


def multiplication_rank_lower_bound(d, e):
    pivots = {}
    for a, b in zip(monomial_labels(d), monomial_labels(e)):
        for total in range(len(a)+len(b)-1):
            row = 0
            for i, index in enumerate(a):
                if 0 <= total-i < len(b):
                    row ^= 1 << (4*index+b[total-i])
            while row:
                pivot = row.bit_length()
                if pivot not in pivots:
                    pivots[pivot] = row
                    break
                row ^= pivots[pivot]
    return len(pivots)


def has_growth_witness(d, e, loss):
    n = len(d)
    for start in range(n):
        for length in range(1, n):
            left = sum((d+d)[start:start+length])-1
            right = 3-sum((e+e)[start:start+length])
            if min(left,right) >= 1 and max(left,right) >= loss+2:
                return True
    return False


def main():
    counts = Counter()
    exceptions = []
    for n in range(2, 9):
        candidates = {orbit(d,e) for d,e in combinations_with_replacement(distributions(n),2)
                      if d != e and all(x+y for x,y in zip(d,e))}
        for d, e in candidates:
            rank = multiplication_rank_lower_bound(d,e)
            witness = has_growth_witness(d,e,8-rank)
            counts[rank,witness] += 1
            if not witness:
                exceptions.append((d,e))
    assert counts == {(7, True):24, (8, True):148, (8, False):2}
    assert set(exceptions) == {orbit((4,0),(3,1)), orbit((3,1,0),(3,0,1))}
    print("PASS: 174 polygon degree pairs; 24 rank-seven lower bounds with a growth witness, 150 full-rank products")


if __name__ == "__main__":
    main()
