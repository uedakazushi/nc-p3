"""Compute the evaluation ranks for a pentagon with five independent gluing ratios."""
from itertools import combinations

from sympy import Matrix, factor, linsolve, ones, symbols, zeros


def main():
    ratios = symbols("q:5", nonzero=True)
    parameters = symbols("a b c")
    edge_entries = (0, *parameters, 0)
    chord_entries = symbols("u:10")
    values = zeros(5)
    for i in range(5):
        neighbor, opposite = (i + 1) % 5, (i + 2) % 5
        values[i, neighbor] = edge_entries[i]
        values[neighbor, i] = -ratios[i]*edge_entries[i]
        values[i, opposite], values[opposite, i] = chord_entries[2*i:2*i+2]

    # Relations evaluated on the frame have zero row and column sums.
    # Their left evaluation at vertex zero vanishes as well.
    constraints = [*values[0, :], *(ones(1, 5)*values), *(values*ones(5, 1))]
    solution, = linsolve(constraints, chord_entries)
    assert not set(solution).intersection(chord_entries)
    values = values.subs(dict(zip(chord_entries, solution)))
    assert not values.free_symbols.intersection(chord_entries)
    assert values[0, :] == zeros(1, 5) and ones(1, 5)*values == zeros(1, 5)
    assert values*ones(5, 1) == zeros(5, 1)

    # Rows 1,2,3 determine the image, since row zero vanishes and rows sum to zero.
    expected = (0, -ratios[1]*ratios[2], -ratios[2], -ratios[3], 1)
    for vertex, certificate in enumerate(expected):
        block = Matrix(values[1:4, vertex]).jacobian(parameters)
        minors = {factor(block.extract(rows, cols).det())
                  for rows in combinations(range(3), 2)
                  for cols in combinations(range(3), 2)}
        if vertex == 0:
            assert minors == {0}
        else:
            assert certificate in minors and block.det().expand() == 0
    print("PASS: exact node ranks <=1, 2, 2, 2, 2 for all nonzero edge ratios")


if __name__ == "__main__":
    main()
