"""Enumerate dimension vectors and both endpoint obstruction bounds with one bounded corner."""

from itertools import product


def possible_vectors():
    survivors = []
    for x, z in product(range(2, 4), repeat=2):
        if z < x:
            continue
        for y in range(2, 3*(x-1)+1):
            d = (2, x, y, z)
            h = [2*d[i] - 4*d[i-1] - 4*d[(i+1)%4] + 6*d[(i+2)%4]
                 for i in range(4)]
            euler = sum(a*b for a, b in zip(d, h))
            lower = []
            for i in range(4):
                left = 3*(d[i-1]-1) - d[i]
                right = 3*(d[(i+1)%4]-1) - d[i]
                bounds = [a*b for a, b in product(range(left+1), range(right+1))
                          if a+b+h[i] >= 0]
                if not bounds:
                    break
                lower.append(min(bounds))
            if len(lower) == 4 and sum(lower) <= euler:
                survivors.append(d)
    return survivors


if __name__ == "__main__":
    assert set(possible_vectors()) == {(2,2,2,2), (2,2,2,3), (2,2,3,2)}
    print("Only the balanced and one exceptional size-three entry survive the first inequalities.")
