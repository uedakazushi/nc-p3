"""Characteristic-zero differential and two transverse marked binary directions."""
from itertools import permutations, product

import numpy as np
import sympy as sp

from verify_binary_component_crossing import (
    character_matrices, reduced, sample, sample_rational, tangent_data,
)
from verify_binary_opposite_determinant import matrix, parity, pivots
from verify_rank_fifteen_opposite_kernels import BOTH, FIRST, PRIME, SECOND


DECIMIC = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
OCTAVIC = [2, 3, 20, 43, 74, 115, 164, 227, 296]
DIRECTION = [3, 1, 4, 1, 5, 9, 2, 6, 5]


def presentation(gram, hankel):
    shifts = [sp.Matrix(6, 4, lambda i, j: int(i == j + k)) for k in range(3)]
    return sp.Matrix.hstack(*shifts, gram.inv()*hankel).T


def marked_presentations():
    gram = sp.Matrix(6, 6, lambda i, j: DECIMIC[i+j])
    hankel = sp.Matrix(6, 4, lambda i, j: OCTAVIC[i+j])
    relation = presentation(gram, hankel)
    x = sp.symbols("x")
    d, linear, cubic = x, x-1, x**3+2

    def column(form, size):
        return sp.Matrix([sp.expand(form).coeff(x, i) for i in range(size)])

    domain = sp.Matrix.hstack(*[column(linear*x**i, 6) for i in range(5)],
                              *[column(cubic*x**i, 6) for i in range(3)])
    image = sp.Matrix.hstack(*[column(d*x**i, 6) for i in range(5)],
                             *[relation[12:, :].T*column(d*x**i, 4) for i in range(3)])
    cols = domain.rref()[1]
    transport = image[:, cols]*domain[:, cols].inv()
    letters = sp.Matrix.hstack(*[column(linear*x**i, 4) for i in range(3)], column(cubic, 4))
    reverse_gram = transport.T*gram*transport
    extra = sp.eye(6)[:, :3].row_join(relation[12, :].T)
    reverse = presentation(reverse_gram, reverse_gram*transport.inv()*extra*letters.inv())
    order = sp.eye(4)[:, [1, 2, 3, 0]].inv()
    change = sp.Matrix.hstack(*[sp.kronecker_product(order[:, j], letters[:, i])
                                for i, j in product(range(4), repeat=2)])
    assert change.T*reverse*reverse_gram*reverse.T*change == relation*gram*relation.T
    return gram, relation, reverse_gram, reverse, change


def deformation(gram, relation, delta_gram, delta_hankel):
    derivative = sp.zeros(16, 6)
    derivative[12:, :] = (gram.inv()*(delta_hankel-delta_gram*relation[12:, :].T)).T
    return (derivative*gram*relation.T + relation*delta_gram*relation.T
            + relation*gram*derivative.T)


def marked_directions():
    gram, relation, reverse_gram, reverse, change = marked_presentations()
    hankel = sp.Matrix(6, 4, lambda i, j: DIRECTION[i+j])
    u = deformation(gram, relation, sp.zeros(6), hankel)
    v = deformation(reverse_gram, reverse, sp.zeros(6), hankel)
    return u, change.T*v*change


def marked_tangent(gram, relation):
    tangent = []
    for j in range(11):
        delta = sp.Matrix(6, 6, lambda i, k: int(i+k == j))
        tangent.append(reduced(deformation(gram, relation, delta, sp.zeros(6, 4))))
    for j in range(9):
        delta = sp.Matrix(6, 4, lambda i, k: int(i+k == j))
        tangent.append(reduced(deformation(gram, relation, sp.zeros(6), delta)))
    w = reduced(relation*gram*relation.T)
    identity = np.eye(4, dtype=int)
    for elementary in np.eye(16, dtype=int).reshape(-1, 4, 4):
        for action in (np.kron(elementary, identity), np.kron(identity, elementary)):
            tangent.append((action.T @ w + w @ action) % PRIME)
    return np.array(tangent)


def transversality():
    gram, relation, reverse_gram, reverse, change = marked_presentations()
    forward = marked_tangent(gram, relation)
    transform = reduced(change)
    reverse = (transform.T @ marked_tangent(reverse_gram, reverse) % PRIME) @ transform % PRIME
    u, v = map(reduced, marked_directions())

    def rank(*arrays):
        return matrix(np.concatenate([a.reshape(-1, 256) for a in arrays])).rank()

    assert (rank(forward), rank(reverse), rank(forward, reverse)) == (43, 43, 44)
    assert rank(forward, v) == rank(reverse, u) == 44


class RationalSystem:
    def __init__(self, tag, kernel):
        self.kernel = kernel
        basis = [sp.Matrix(v) for v in character_matrices(*tag)]
        self.vectors = sp.Matrix.vstack(*[m.reshape(1, 256) for m in basis])
        self.indices = list(zip(*np.triu_indices(10, int(tag[0] != tag[1]))))
        self.jacobian = sp.Matrix.hstack(*[self.project(m) for m in basis])
        residue = reduced(self.jacobian)
        self.cols = pivots(residue)
        self.rows = pivots(residue[:, self.cols].T)
        self.inverse = self.jacobian[self.rows, self.cols].inv()

    def project(self, value):
        normal = self.kernel.T*value*self.kernel
        return sp.Matrix([normal[i, j] for i, j in self.indices])

    def solve(self, value):
        target = self.project(value)
        coeff = self.inverse*target[self.rows, :]
        assert self.jacobian[:, self.cols]*coeff == target
        return (coeff.T*self.vectors[self.cols, :]).reshape(16, 16)

    def lift(self, residue):
        for coeff in self.jacobian.nullspace():
            tangent = (coeff.T*self.vectors).reshape(16, 16)
            if np.array_equal(reduced(tangent), residue):
                return tangent
        raise ValueError("No rational tangent with the specified reduction")


def exact_differential():
    w, kernel, inverse = sample_rational(OCTAVIC, 1)
    _, integral, n, j = sample(OCTAVIC, 1)
    _, mixed, _ = tangent_data(integral, n, j)
    systems = {tag: RationalSystem(tag, kernel) for tag in (FIRST, SECOND, BOTH)}
    a, b = (systems[tag].lift(mixed[tag]) for tag in (FIRST, SECOND))

    def cross(x, y):
        return x*inverse*y+y*inverse*x

    def opposite(m):
        return sp.Matrix(16, 16, lambda i, j: m[4*(i//4)+j//4, 4*(i%4)+j%4])

    sy, sk = sp.Matrix(parity(1)), sp.Matrix(parity(-1))

    def cut(even, first, second, both):
        top = (sy.T*opposite(even)*sy).row_join(sy.T*opposite(first)*sk)
        bottom = (sk.T*opposite(second)*sy).row_join(sk.T*opposite(both)*sk)
        return top.col_join(bottom)

    z = systems[BOTH].solve(cross(a, b))
    base = cut(w, a, b, z)
    left, right = [sp.Matrix.hstack(*m.nullspace()) for m in (base.T, base)]
    assert left.cols == right.cols == 2
    for direction in marked_directions():
        first = systems[FIRST].solve(cross(direction, a))
        second = systems[SECOND].solve(cross(direction, b))
        source = cross(direction, z)+cross(first, b)+cross(a, second)
        for x, y, v in permutations((direction, a, b)):
            source -= x*inverse*y*inverse*v
        derivative = cut(direction, first, second, systems[BOTH].solve(source))
        assert left.T*derivative*right == sp.zeros(2)


if __name__ == "__main__":
    transversality()
    exact_differential()
    print("Transverse binary markings, rational scaled rank fourteen and both zero axial differentials verified.")
