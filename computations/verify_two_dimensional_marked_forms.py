"""Dimensions of all simple-marked tensor spaces and the two nontrivial cycle conditions."""
from itertools import combinations, product

import sympy as sp


def subsets(positions):
    return [set(part) for n in range(len(positions)+1) for part in combinations(sorted(positions), n)]


def term_dimension(deficient, killed):
    remaining = set(range(4)) - killed
    adjacent = {dart: set() for i in remaining for dart in (2*i, 2*i+1)}

    def join(a, b):
        adjacent[a].add(b)
        adjacent[b].add(a)

    for i in remaining:
        if i in deficient:
            join(2*i, 2*i+1)
        successor = (i+1) % 4
        if successor in remaining:
            join(2*i+1, 2*successor)
    pending = set(adjacent)
    dimension = 1
    while pending:
        chain = {pending.pop()}
        frontier = set(chain)
        while frontier:
            neighbors = set().union(*(adjacent[dart] for dart in frontier)) & pending
            chain.update(neighbors)
            pending.difference_update(neighbors)
            frontier = neighbors
        dimension *= len(chain) + 1
    return dimension


def twisted_symmetric_dimension(transform):
    # A symmetric binary tensor entry depends only on the number of ones.
    block = sp.kronecker_product(transform, transform.inv()) - sp.eye(4)
    rows = []
    for count in range(7):
        for row in range(4):
            coefficients = [0]*9
            for a, b in product(range(2), repeat=2):
                coefficients[count+a+b] += block[row, 2*a+b]
            rows.append(coefficients)
    return 9 - sp.Matrix(rows).rank()


def main():
    for deficient in subsets(range(4)):
        total = sum(term_dimension(deficient, killed) for killed in subsets(deficient))
        assert total == (88 if len(deficient) == 4 else 81)
    for matrix, expected in [(sp.eye(2), 9), (sp.diag(1, -1), 2), (sp.Matrix([[1,1],[0,1]]), 2)]:
        dimension = twisted_symmetric_dimension(matrix)
        assert dimension == expected
        assert 88 - 9 + dimension == (88 if expected == 9 else 81)
    print("All sixteen arrow-rank masks and the three holonomy types verified.")


if __name__ == "__main__":
    main()
