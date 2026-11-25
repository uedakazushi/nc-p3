"""Exhaust all invariant Fermat sextic polarizations using exact good reduction."""
from collections import Counter
from functools import cache
from itertools import product

from flint import nmod_mat, nmod_poly


PRIME = 101
BRANCH = (1, 10, 100, 91)
GROUP = tuple(product(range(4), repeat=2))
LABELS = ((1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2))
T = nmod_poly([0, 1], PRIME)
RAMIFICATION = T**4 - 1


def basis(label):
    m, n = label
    divisor = (6-m-n, m, n, 0)
    sections = []
    for character in range(4):
        poles = [(a+character)//4 for a in divisor]
        for degree in range(sum(poles)-character+1):
            sections.append((character, degree, poles))
    assert len(sections) == 4
    return sections


@cache
def relations(a, b):
    products = []
    for (i, k, x), (j, l, y) in product(basis(a), basis(b)):
        power, character = divmod(i+j, 4)
        products.append((character, power, k+l, [v+w for v, w in zip(x, y)]))
    rows = []
    for character in range(4):
        group = [p for p in products if p[0] == character]
        bounds = [max(p[3][u] for p in group) for u in range(4)]
        coefficients = []
        for e, power, degree, poles in products:
            value = nmod_poly([], PRIME)
            if e == character:
                value = T**degree * RAMIFICATION**power
                for root, bound, pole in zip(BRANCH, bounds, poles):
                    value *= (T-root)**(bound-pole)
            coefficients.append(value)
        rows.extend([[int(v[h]) for v in coefficients]
                     for h in range(1+max(v.degree() for v in coefficients))])
    nullspace, dim = nmod_mat(rows, PRIME).nullspace()
    assert dim == 6
    return [[int(nullspace[i, j]) for i in range(16)] for j in range(dim)]


@cache
def cubic_overlap(a, b, c):
    words = list(product(range(4), repeat=3))
    rows = [[u[4*i+j] if k == v else 0 for i, j, k in words]
            for u in relations(a, b) for v in range(4)]
    rows += [[u[4*j+k] if i == v else 0 for i, j, k in words]
             for u in relations(b, c) for v in range(4)]
    return 48-nmod_mat(rows, PRIME).rank()


def main():
    assert all(pow(root, 4, PRIME) == 1 for root in BRANCH)
    fixed = ((0, 0), (1, 0), (0, 1), (3, 3))
    effective = {tuple((i+j) % 4 for i, j in zip(a, b)) for a, b in product(fixed, repeat=2)}
    assert set(GROUP) - effective == set(LABELS)
    counts = Counter()
    for start, step, curvature in product(LABELS, GROUP, product((0, 2), repeat=2)):
        seq = tuple(tuple((v+i*d+i*(i-1)//2*k) % 4 for v, d, k in zip(start, step, curvature))
                    for i in range(4))
        if all(seq[i] in LABELS and seq[i] != seq[i-1] for i in range(4)):
            counts[tuple(cubic_overlap(seq[i], seq[(i+1) % 4], seq[(i+2) % 4])
                         for i in range(4))] += 1
    assert counts == {(0,)*4: 6, (1,)*4: 16, (2,)*4: 8}
    print("PASS: all 30 ordered polarization windows excluded; cubic overlap bounds 0, 1, 2")


if __name__ == "__main__":
    main()
