"""Certify four unique normal triple-point octics on the rational nonsymmetric chart."""

from itertools import combinations, combinations_with_replacement, permutations
from math import prod

from flint import nmod_mat, nmod_mpoly_ctx

from verify_configuration_intersections import (
    PAIRS, change_of_square_coordinates, exchange, modular_points, partner, rational_example,
)


def exponents(degree):
    return [(i, j, degree-i-j) for i in range(degree+1) for j in range(degree+1-i)]


def plane_section(change, prime):
    ctx = nmod_mpoly_ctx.get(("x", "y", "z"), prime)
    x, y, z = ctx.gens()
    letters = [x, y, z, x+2*y+5*z]
    entries = [sum(int(change[k, j])*letters[a]*letters[b] for k, (a, b) in enumerate(PAIRS))
               for j in range(10)]
    result = ctx.constant(0)
    for perm in permutations(range(4)):
        sign = (-1)**sum(perm[i] > perm[j] for i, j in combinations(range(4), 2))
        result += sign*prod(entries[PAIRS.index(tuple(sorted((i, perm[i]))))] for i in range(4))
    smooth_plane_form(result, prime)


def smooth_plane_form(result, prime):
    rows = []
    for i in range(3):
        derivative = result.derivative(i).to_dict()
        for multiplier in exponents(12):
            rows.append([int(derivative.get(tuple(a-b for a, b in zip(monomial, multiplier)), 0))
                         for monomial in exponents(19)])
    assert nmod_mat(rows, prime).rank() == 210


def triple_evaluation(points, prime):
    words = list(combinations_with_replacement(range(4), 8))
    rows = []
    for point in points:
        for i, j in PAIRS:
            row = []
            for word in words:
                coefficient = word.count(i)*(word.count(j)-int(i == j))
                if coefficient:
                    remainder = list(word)
                    remainder.remove(i)
                    remainder.remove(j)
                    coefficient *= prod(int(point[k]) for k in remainder)
                row.append(coefficient % prime)
            rows.append(row)
    assert nmod_mat(rows, prime).rank() == 164


def main():
    tensor, _, _ = rational_example()
    prime = 2657
    opposite = exchange(tensor)
    a = modular_points(tensor, prime)
    b = partner(a, tensor, False, prime)
    c = partner(b, opposite, True, prime)
    d = partner(c, tensor.T, False, prime)
    for first, third in [(a, c), (b, d), (c, a), (d, b)]:
        triple_evaluation(first, prime)
        plane_section(change_of_square_coordinates(first, third, prime), prime)
    print("All four jet ranks are 164 and each determinantal octic has a smooth plane section.")


if __name__ == "__main__":
    main()
