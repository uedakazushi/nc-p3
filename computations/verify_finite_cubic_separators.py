"""Check the four-step cubic separator in an opposite-rank-ten, finite example."""
from itertools import product

from sympy import I, Matrix, eye, simplify, sqrt, symbols, zeros

from verify_twenty_point_companion import gram_half


def main():
    gram = gram_half(1).T * gram_half(-1)
    t, x, y, z = variables = symbols("t x y z")
    matrix = Matrix(4, 4, lambda j, l: sum(variables[k]*gram[j, 4*k+l] for k in range(4)))
    assert matrix[0, :] == zeros(1, 4) and matrix[:, 0] == zeros(4, 1)
    compressed = matrix[1:, 1:]
    assert compressed.det().expand() == (t*(t*t+x*x+y*y+z*z)).expand()
    assert Matrix.hstack(*(compressed.diff(v).reshape(9, 1) for v in variables)).rank() == 4
    points = [(1, 0, 0, 0)]
    for axis in range(3):
        points.append((0, *(int(i == axis) for i in range(3))))
        for signs in product((-1, 1), repeat=2):
            nonzero = iter(signs)
            points.append((1, *(0 if i == axis else next(nonzero)*I/sqrt(2) for i in range(3))))
    points.extend((0, 1, a, b) for a, b in product((-1, 1), repeat=2))
    assert len(points) == 20
    assert [simplify(compressed.det().subs(dict(zip(variables, p)))) for p in points] == [1]+[0]*19
    assert compressed.subs(dict(zip(variables, points[0]))) == eye(3)
    print("PASS: four independent coefficients, full-rank matched value and all nineteen separator zeros")


if __name__ == "__main__":
    main()
