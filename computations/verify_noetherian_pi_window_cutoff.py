"""Check the alternating matrix-unit identity and the uniform central degrees."""
from itertools import permutations, product
from math import lcm


def signed_matrix_unit_value(xs, separators, first_row, last_column):
    size = len(xs)
    value = 0
    for order in permutations(range(size)):
        path = [xs[order[0]]]
        for separator, position in zip(separators, order[1:]):
            path.extend([separator, xs[position]])
        if path[0][0] != first_row or path[-1][1] != last_column:
            continue
        if all(left[1] == right[0] for left, right in zip(path, path[1:])):
            inversions = sum(i > j for idx, i in enumerate(order) for j in order[idx + 1:])
            value += (-1)**inversions
    return value


def main():
    xs = list(product(range(2), repeat=2))
    separators = [(left[1], right[0]) for left, right in zip(xs, xs[1:])]
    component = signed_matrix_unit_value(xs, separators, 0, 1)
    assert component == 1
    diagonal = [3, -7]
    insertions = sum(diagonal[row] * component for row, _ in xs)
    assert insertions == 2 * sum(diagonal) * component
    for bound in range(1, 17):
        degree = 4 * lcm(*range(1, bound + 1)) * lcm(*range(1, bound*bound + 1))
        for i, j in product(range(1, bound + 1), range(1, bound*bound + 1)):
            assert degree % (4*i*j) == 0
    print('Alternating matrix-unit component, left-trace multiplier, and parameter divisibilities verified.')


if __name__ == '__main__':
    main()
