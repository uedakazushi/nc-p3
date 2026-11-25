"""Check the linkage numbers, Picard identity and cubic overlap of one alternating pair."""
from itertools import product
from math import comb

from flint import fmpq_mat
from sympy import Poly, symbols


def overlap(left, right):
    words = list(product(range(4), repeat=3))
    rows = [[rel[4*i+j] if k == a else 0 for i, j, k in words]
            for rel in left for a in range(4)]
    rows += [[rel[4*j+k] if i == a else 0 for i, j, k in words]
             for rel in right for a in range(4)]
    return len(rows) - fmpq_mat(rows).rank()


def mixed_row(quadric, variables):
    row = [0] * 16
    for powers, coefficient in Poly(quadric, *variables).terms():
        i, j = [k for k, power in enumerate(powers) for _ in range(power)]
        if i == 4 or j == 3:
            i, j = j, i
        assert i < 4 and j != 3
        row[4*i + (3 if j == 4 else j)] += int(coefficient)
    return row


def main():
    t = symbols("t")
    assert ((1-t)**3 * (1+t)).expand() == 1-2*t+2*t**3-t**4
    assert -10 - (-6) == 1-5
    assert all((a-1)*(6-a) != 5 for a in range(8))
    for n in range(25):
        assert sum(comb(h+2, 2) for h in range(n+1)) == comb(n+3, 3)
    h, k, z, x, y = variables = symbols("h k z x y")
    r, s, u, v = x+y+h, x+y-h, y-x+k, y-x-k
    net = (r*r+2*s*s+3*u*u-6*v*v, 2*r*r-3*s*s+5*u*u-4*v*v,
           z*z-h*(x+2*y)-k*(3*x+5*y)-h*h-3*k*k-2*h*k)
    relations = []
    for a, b in ((0, 1), (0, 2), (1, 2)):
        row = [0] * 16
        row[4*a+b], row[4*b+a] = 1, -1
        relations.append(row)
    relations.extend(mixed_row(q, variables) for q in net)
    opposite = [[row[4*j+i] for i, j in product(range(4), repeat=2)]
                for row in relations]
    assert overlap(relations, opposite) == overlap(opposite, relations) == 3
    print("PASS: linkage genus, polarization identity, polynomial-module counts,"
          " and rejection of a three-overlap alternating pair")


if __name__ == "__main__":
    main()
