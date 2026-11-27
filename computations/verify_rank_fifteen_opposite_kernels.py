"""Fourth-order Schur calculation for the two bilinear kernels of the residual cut."""
import numpy as np

from verify_binary_component_crossing import (
    PRIME, character_matrices, mod_kernel, reduced, sample, tangent_data,
)
from verify_binary_opposite_determinant import matrix, multiply, opposite, parity, pivots


PLUS, FIRST, SECOND, BOTH = (1, 1), (1, -1), (-1, 1), (-1, -1)


def modular_inverse(array):
    return np.array(matrix(array).inv().tolist(), dtype=int)


class CharacterSystem:
    def __init__(self, kernel, character):
        self.kernel = kernel
        self.basis = character_matrices(*character)
        self.index = np.triu_indices(10, int(character[0] != character[1]))
        normal = multiply(kernel.T, self.basis, kernel)
        self.jacobian = normal[:, self.index[0], self.index[1]].T
        self.cols = pivots(self.jacobian)
        self.rows = pivots(self.jacobian[:, self.cols].T)
        self.inverse = modular_inverse(self.jacobian[np.ix_(self.rows, self.cols)])
        self.cokernel = mod_kernel(self.jacobian.T)

    def target(self, source):
        return multiply(self.kernel.T, source, self.kernel)[self.index]

    def obstruction(self, source):
        return multiply(self.cokernel.T, self.target(source))

    def solve(self, source):
        target = self.target(source)
        values = multiply(self.inverse, target[self.rows])
        assert np.array_equal(multiply(self.jacobian[:, self.cols], values), target)
        return (values @ self.basis[self.cols].reshape(-1, 256)).reshape(16, 16) % PRIME


def schur_source(degree, layers, inverse):
    """Nonlinear terms moved to the right hand side of the Schur equation."""
    result = {}

    def extend(remaining, term, character, sign):
        if remaining == 0:
            result[character] = (result.get(character, 0) + sign * term) % PRIME
            return
        for index, layer in enumerate(layers[:remaining], start=1):
            for tag, value in layer.items():
                new_tag = (character[0] * tag[0], character[1] * tag[1])
                product = value if term is None else multiply(term, inverse, value)
                extend(remaining - index, product, new_tag, -sign)

    extend(degree, None, PLUS, 1)
    return result


def scaled_cut(even, first, second, both):
    sy, sk = parity(1), parity(-1)
    return np.block([[multiply(sy.T, opposite(even), sy), multiply(sy.T, opposite(first), sk)],
                     [multiply(sk.T, opposite(second), sy), multiply(sk.T, opposite(both), sk)]])


def main():
    w, integral, kernel, inverse = sample()
    symmetric, mixed, forms = tangent_data(integral, kernel, inverse)
    systems = {tag: CharacterSystem(kernel, tag) for tag in (PLUS, FIRST, SECOND, BOTH)}
    linear = {tag: mixed[tag] for tag in (FIRST, SECOND)}
    source = schur_source(2, [linear], inverse)
    even = systems[PLUS].solve(source[PLUS])
    both = systems[BOTH].solve(source[BOTH])
    leading = scaled_cut(reduced(w), linear[FIRST], linear[SECOND], both)
    assert matrix(leading).rank() == 15
    left, right = mod_kernel(leading.T)[:, 0], mod_kernel(leading)[:, 0]
    sy, sk = parity(1), parity(-1)
    assert not left[10:].any()
    cone = (sy @ left[:10] % PRIME).reshape(4, 4)
    assert not cone[3, :].any() and not cone[:, 3].any()
    assert matrix(cone[:3, :3]).rank() == 3
    assert matrix((sk @ right[10:] % PRIME).reshape(4, 4)).rank() == 4

    def cubic_obstruction(value):
        cubic = schur_source(3, [linear, {PLUS: value, BOTH: both}], inverse)
        return int(systems[FIRST].obstruction(cubic[FIRST])[0])

    obstruction = cubic_obstruction(even)
    if obstruction:
        for tangent in symmetric:
            derivative = (cubic_obstruction((even + tangent) % PRIME) - obstruction) % PRIME
            if derivative:
                even = (even - obstruction * pow(derivative, -1, PRIME) * tangent) % PRIME
                break
    assert cubic_obstruction(even) == 0

    def quadratic_cut(value):
        layers = [linear, {PLUS: value, BOTH: both}]
        cubic = schur_source(3, layers, inverse)
        layers.append({tag: systems[tag].solve(cubic[tag]) for tag in (FIRST, SECOND)})
        quartic = schur_source(4, layers, inverse)
        return scaled_cut(value, layers[2][FIRST], layers[2][SECOND], systems[BOTH].solve(quartic[BOTH]))

    correction = quadratic_cut(even)
    defect = int(multiply(left, correction, right))
    tangents = mod_kernel(forms).T @ symmetric.reshape(len(symmetric), 256) % PRIME
    for tangent in tangents.reshape(-1, 16, 16):
        trial = quadratic_cut((even + tangent) % PRIME)
        derivative = (int(multiply(left, trial, right)) - defect) % PRIME
        if derivative:
            even = (even - defect * pow(derivative, -1, PRIME) * tangent) % PRIME
            correction = quadratic_cut(even)
            break
    assert multiply(left, correction, right) == 0

    rows, cols = pivots(leading.T), pivots(leading)
    target = -multiply(left, correction) % PRIME
    lift = np.zeros(16, dtype=int)
    lift[rows] = multiply(target[cols], modular_inverse(leading[np.ix_(rows, cols)]))
    assert np.array_equal(multiply(lift, leading), target)
    symmetric_lift = (sy @ lift[:10] % PRIME).reshape(4, 4)
    skew_lift = (sk @ lift[10:] % PRIME).reshape(4, 4)
    product = multiply(skew_lift[3, :3], modular_inverse(cone[:3, :3]), skew_lift[:3, 3])
    coefficient = int((symmetric_lift[3, 3] - product) % PRIME)
    assert coefficient == 22
    print("Right skew rank four and left determinant Schur coefficient 22 modulo 1009 verified.")


if __name__ == "__main__":
    main()
