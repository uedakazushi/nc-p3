"""Test every cubic point contraction on the symmetric twenty-point example."""
from itertools import combinations, combinations_with_replacement, product

from sympy import Matrix, groebner, ilcm, symbols

from verify_finite_square_chart import frame


def main():
    quartics = list(combinations_with_replacement(range(4), 4))
    equations = []
    for a, (b, c) in product(range(4), combinations_with_replacement(range(4), 2)):
        row = [0]*len(quartics)
        row[quartics.index(tuple(sorted((a, a, b, c))))] += 1
        for i, j in product(range(4), repeat=2):
            row[quartics.index(tuple(sorted((i, j, b, c))))] -= 1
        equations.append(row)
    socle, = Matrix(equations).nullspace()
    denominator = ilcm(*(v.q for v in socle))
    tensor = dict(zip(quartics, denominator * socle))
    squares = list(combinations_with_replacement(range(3), 2))
    cubes = list(combinations_with_replacement(range(3), 3))
    for pair in combinations(squares, 2):
        assert sum(all(divisor in pair for divisor in combinations(word, 2))
                   for word in cubes) <= 2
    points, _ = frame()
    for point in points.tolist():
        def entry(i, j, k):
            return sum(v*tensor[tuple(sorted((a, i, j, k)))] for a, v in enumerate(point))
        slices = [Matrix(4, 4, lambda j, k: entry(i, j, k)) for i in range(4)]
        span = Matrix.vstack(*(m.reshape(1, 16) for m in slices))
        _, basis = span.T.rref()
        assert len(basis) == 3
        coords = symbols("x y z")
        pencil = sum((v*slices[i][basis, basis] for v, i in zip(coords, basis)), Matrix.zeros(3))
        minors = [pencil[ij, kl].det() for ij in combinations(range(3), 2)
                  for kl in combinations(range(3), 2)]
        for i, x in enumerate(coords):
            chart = {x: 1, **{y: 0 for y in coords[:i]}}
            values = [f.subs(chart) for f in minors]
            if i < 2:
                assert list(groebner(values, *coords[i+1:], modulus=101)) == [1]
            else:
                assert any(v % 101 for v in values)
    print("PASS: all twenty cubic contractions have three supports of dimension 3"
          " and rank-one-free ternary nets")


if __name__ == "__main__":
    main()
