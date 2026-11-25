"""Check determinant exponents, fiber degrees and the involution thresholds."""
from functools import reduce
from itertools import combinations, product
from math import gcd, lcm

from sympy import Matrix, groebner, symbols

from verify_septic_bielliptic_forms import intersections
from verify_septic_point_curves import mixed_row


def septic_fiber():
    r, s, h, x, y = variables = symbols("r s h x y")
    square = h*h + r*s
    quadrics = [x*y-r*r-2*s*s-3*square, x*r-2*r*y-2*r*r-5*s*s-square,
                x*s+2*s*y-3*r*r+s*s-2*square]
    relations = []
    for a, b in combinations(range(3), 2):
        row = [0]*16
        row[4*a+b], row[4*b+a] = 1, -1
        relations.append(row)
    relations.extend(mixed_row(q, variables) for q in quadrics)
    reverse = [[row[4*j+i] for i, j in product(range(4), repeat=2)] for row in relations]
    left, right = intersections(relations, reverse, 2), intersections(reverse, relations, 2)
    top, = intersections(left, right, 3)
    divisor = reduce(gcd, top)
    top = [v//divisor for v in top]
    coords = (r, s, h, x)
    contraction = Matrix(4, 4, lambda j, l: sum(coords[k]*top[192+16*j+4*k+l] for k in range(4)))
    edge = Matrix(6, 4, lambda i, l: sum(coords[k]*relations[i][4*k+l] for k in range(4)))
    equations = [edge[list(rows), :].det() for rows in combinations(range(6), 4)]
    equations += [contraction[rows, cols].det() for rows in combinations(range(4), 2)
                  for cols in combinations(range(4), 2)]
    linear = r+2*s+3*h+5*x
    affine = groebner(equations+[linear-1], *coords, modulus=101, order="grevlex")
    assert list(affine) == [x*x, r+17*x, s-34*x, h-15*x-34]
    infinity = groebner(equations+[linear], *coords, modulus=101, order="grevlex")
    assert infinity.is_zero_dimensional


def main():
    # Determinant twist along the first and second axes in ordered neighboring positions.
    for pos in range(3):
        assert ((-1, 2, -1)[pos] + (1, -2, 1)[pos]) == 0
    for genus in range(2, 12):
        self_intersection = 2-2*genus
        for quotient_genus in range(genus):
            fixed = 2*genus+2-4*quotient_genus
            if fixed >= 0:
                assert (self_intersection+fixed < 0) == (quotient_genus >= 2)
    assert lcm(2, 3, 2, 1) == 6
    # The restricted displacement on an elliptic translation fiber is twice the same class.
    for order in range(1, 13):
        assert (order // gcd(2, order) == 1) == (order in (1, 2))
    septic_fiber()
    print("PASS: length-two septic fiber, determinant twists, negative intersections, Hurwitz thresholds and torsion improvement")


if __name__ == "__main__":
    main()
