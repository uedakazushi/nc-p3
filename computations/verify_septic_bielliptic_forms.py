"""Smooth characteristic-zero examples and exact overlap ranks in both parities."""
from itertools import combinations, product

from flint import fmpz_mat
from sympy import Matrix, groebner, symbols

from verify_septic_point_curves import mixed_row


def smooth_projective(equations, coordinates):
    jacobian = Matrix(equations).jacobian(coordinates)
    minors = [jacobian[:, cols].det() for cols in combinations(range(4), len(equations))]
    for i, coord in enumerate(coordinates):
        chart = {coord: 1, **{v: 0 for v in coordinates[:i]}}
        equations_chart = [f.subs(chart) for f in equations + minors]
        if i == len(coordinates) - 1:
            assert any(int(f) % 101 for f in equations_chart)
        else:
            assert list(groebner(equations_chart, *coordinates[i+1:],
                                 order="grevlex", modulus=101)) == [1]


def intersections(left, right, degree):
    size = 4 ** (degree + 1)
    rows = [[v[j // 4] if j % 4 == a else 0 for j in range(size)]
            for v in left for a in range(4)]
    cut = len(rows)
    rows += [[v[j % (size // 4)] if j // (size // 4) == a else 0 for j in range(size)]
             for a in range(4) for v in right]
    dependencies, rank = fmpz_mat(rows).transpose().nullspace()
    return [[sum(int(dependencies[i, k]) * c for i, c in enumerate(column))
             for column in zip(*rows[:cut])] for k in range(rank)]


def flatten_ranks(tensor):
    indices = list(product(range(4), repeat=4))
    for i in range(4):
        rotated = sorted(zip(indices, tensor), key=lambda pair: pair[0][i:] + pair[0][:i])
        entries = [value for _, value in rotated]
        for rows, cols, rank in ((4, 64, 4), (16, 16, 6)):
            assert fmpz_mat(rows, cols, entries).rank() == rank


def verify(parameter):
    r, s, h, x, y = variables = symbols("r s h x y")
    square = h*h + parameter*r*s
    quadrics = [x*y-r*r-2*s*s-3*square, x*r-2*r*y-2*r*r-5*s*s-square,
                x*s+2*s*y-3*r*r+s*s-2*square]
    base = [(3*q-c*quadrics[0]).expand() for q, c in zip(quadrics[1:], (1, 2))]
    smooth_projective(base, (r, s, x, y))
    smooth_projective([q.subs(h, 0) for q in quadrics], (r, s, x, y))
    relations = []
    for a, b in combinations(range(3), 2):
        row = [0]*16
        row[4*a+b], row[4*b+a] = 1, -1
        relations.append(row)
    relations.extend(mixed_row(q, variables) for q in quadrics)
    reverse = [[row[4*j+i] for i, j in product(range(4), repeat=2)] for row in relations]
    first = intersections(relations, reverse, 2)
    second = intersections(reverse, relations, 2)
    assert len(first) == len(second) == 4
    for left, right in ((first, second), (second, first)):
        top = intersections(left, right, 3)
        assert len(top) == 1
        flatten_ranks(top[0])


if __name__ == "__main__":
    for parameter in (0, 1):
        verify(parameter)
    print("PASS: smooth base and branch modulo 101, exact rational ranks 4,6,4,1 in both parities")
