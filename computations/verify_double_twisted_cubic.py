"""Eliminate the cyclic tensor of a canonically embedded nonsplit genus-three ribbon."""
from collections import defaultdict
from itertools import product
from math import comb

from sympy import I, Matrix, Rational, expand, groebner, symbols

a, b, c = symbols("a b c")


def first_edge_equations(word):
    """Multiply four translated brackets, differentiate on the ribbon, restrict to y=x+1."""
    terms = {word: 1}
    for i in range(4):
        result = defaultdict(int)
        for powers, coefficient in terms.items():
            result[powers] -= coefficient
            for j, sign in ((i, -1), ((i+1) % 4, 1)):
                new = list(powers)
                new[j] += 1
                result[tuple(new)] += sign*coefficient
        terms = result

    equations = defaultdict(int)
    for powers, coefficient in terms.items():
        for slot in (0, 1):
            k = powers[slot]
            for drop, scalar in enumerate((2*a, 2*b+2*k-3, 2*c), start=1):
                if drop > k:
                    continue
                new = list(powers)
                new[slot] -= drop
                for n in range(new[1]+1):
                    equations[(new[0]+n, new[2], new[3])] += (
                        coefficient*scalar*comb(new[1], n))
    return equations


def rotation_matrix(eigenvalue):
    columns = []
    seen = set()
    for word in product(range(2), repeat=4):
        if word in seen:
            continue
        orbit = [word]
        rotated = word[1:] + word[:1]
        while rotated != word:
            orbit.append(rotated)
            rotated = rotated[1:] + rotated[:1]
        seen.update(orbit)
        if eigenvalue**len(orbit) != 1:
            continue
        column = defaultdict(int)
        for j, member in enumerate(orbit):
            for exponent, value in first_edge_equations(member).items():
                column[exponent] += eigenvalue**j*value
        columns.append(column)
    keys = sorted(set().union(*(set(column) for column in columns)))
    rows = [[expand(column[key]) for column in columns] for key in keys]
    return Matrix([row for row in rows if any(row)])


def main():
    t = symbols("t")
    coordinates = [t**i for i in range(4)]
    correction = [0, a, a*t+b+Rational(1, 2), a*t**2+(b+Rational(3, 2))*t+c]
    conormals = [expand(coordinates[i]*correction[j+1]+correction[i]*coordinates[j+1]
                        -coordinates[i+1]*correction[j]-correction[i+1]*coordinates[j])
                 for i, j in ((0, 1), (0, 2), (1, 2))]
    block = Matrix([[p.coeff(t, i) for i in range(3)] for p in conormals])
    assert expand(block.det() - (Rational(1, 4)-b**2+a*c)) == 0

    for root, size in ((1, 6), (-1, 4), (I, 3)):
        matrix = rotation_matrix(root)
        assert matrix.cols == size
        coefficients = symbols("f:"+str(size))
        equations = matrix*Matrix(coefficients)
        for first in range(size):
            fixed = {u: int(i == first) for i, u in enumerate(coefficients[:first+1])}
            chart = [expand(eq.subs(fixed)) for eq in equations]
            basis = groebner([eq for eq in chart if eq != 0],
                             *coefficients[first+1:], a, b, c,
                             order="grevlex", extension=I, method="f5b")
            if root == 1 and first == 0:
                assert basis.reduce(a**2)[1] == 0
                assert basis.reduce((b+Rational(1, 2))**2)[1] == 0
            else:
                assert [p.as_expr() for p in basis.polys] == [1]
    print("PASS: conormal determinant and all projective cyclic charts over Q(i); conjugation covers -i")


if __name__ == "__main__":
    main()
