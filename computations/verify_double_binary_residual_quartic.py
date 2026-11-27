"""The transverse quartic of the residual determinant at a double binary marking."""
import numpy as np

from verify_binary_component_crossing import mod_kernel, reduced, sample, tangent_data
from verify_binary_opposite_determinant import matrix, multiply, pivots
from verify_double_binary_transversal import (
    OCTAVIC, exact_differential, marked_directions, transversality,
)
from verify_rank_fifteen_opposite_kernels import (
    BOTH, FIRST, PLUS, PRIME, SECOND, CharacterSystem, modular_inverse, scaled_cut, schur_source,
)


def main():
    transversality()
    exact_differential()
    w, integral, n, inverse = sample(OCTAVIC, 1)
    symmetric, mixed, forms = tangent_data(integral, n, inverse)
    systems = {tag: CharacterSystem(n, tag) for tag in (PLUS, FIRST, SECOND, BOTH)}
    linear = {tag: mixed[tag] for tag in (FIRST, SECOND)}
    source = schur_source(2, [linear], inverse)
    both = systems[BOTH].solve(source[BOTH])
    leading = scaled_cut(reduced(w), linear[FIRST], linear[SECOND], both)
    assert matrix(leading).rank() == 14
    left, right = mod_kernel(leading.T), mod_kernel(leading)
    rows, cols = pivots(leading.T), pivots(leading)
    inverse_cut = np.zeros((16, 16), dtype=int)
    inverse_cut[np.ix_(cols, rows)] = modular_inverse(leading[np.ix_(rows, cols)])

    def quadratic(direction, a, b, check_tangents=False):
        first = {PLUS: direction, FIRST: a*mixed[FIRST] % PRIME, SECOND: b*mixed[SECOND] % PRIME}
        row_scale = np.diag([1]*10+[pow(b, -1, PRIME)]*6)
        col_scale = np.diag([1]*10+[pow(a, -1, PRIME)]*6)

        def cut(even, first, second, both):
            return multiply(row_scale, scaled_cut(even, first, second, both), col_scale)

        source = schur_source(2, [first], inverse)
        second = {tag: systems[tag].solve(value) for tag, value in source.items()}

        def cubic(even):
            layer = {**second, PLUS: even}
            return schur_source(3, [first, layer], inverse)

        def obstruction(even):
            return int(systems[FIRST].obstruction(cubic(even)[FIRST])[0])

        even = second[PLUS]
        defect = obstruction(even)
        if defect:
            for tangent in symmetric:
                slope = (obstruction((even+tangent) % PRIME)-defect) % PRIME
                if slope:
                    even = (even-defect*pow(slope, -1, PRIME)*tangent) % PRIME
                    break

        def projected(even):
            layer = {**second, PLUS: even}
            third = {tag: systems[tag].solve(value) for tag, value in cubic(even).items()}
            fourth = schur_source(4, [first, layer, third], inverse)
            linear_cut = cut(direction, layer[FIRST], layer[SECOND], third[BOTH])
            assert not multiply(left.T, linear_cut, right).any()
            quadratic_cut = cut(even, third[FIRST], third[SECOND], systems[BOTH].solve(fourth[BOTH]))
            return multiply(left.T, (quadratic_cut-multiply(linear_cut, inverse_cut, linear_cut)) % PRIME, right)

        value = projected(even)
        if check_tangents:
            tangents = mod_kernel(forms).T @ symmetric.reshape(len(symmetric), 256) % PRIME
            for tangent in tangents.reshape(-1, 16, 16):
                assert np.array_equal(projected((even+tangent) % PRIME), value)
        return value

    u, v = map(reduced, marked_directions())
    zero = np.zeros((16, 16), dtype=int)
    offset = quadratic(zero, 1, 1, check_tangents=True)
    assert offset.tolist() == [[0, 546], [570, 22]]
    for x, y, a, b in ((1,0,1,1), (0,1,1,1), (1,1,1,1), (0,0,2,1), (0,0,1,2), (3,5,2,4)):
        value = quadratic((x*u+y*v) % PRIME, a, b)
        assert value[0, 0] == 0
        assert value[0, 1] == (664*y*y+546*b*b) % PRIME
        assert value[1, 0] == (691*x*x+570*a*a) % PRIME
    print("The residual tangent quartic reduces, up to a unit, to (664v²+546b²)(691u²+570a²).")


if __name__ == "__main__":
    main()
