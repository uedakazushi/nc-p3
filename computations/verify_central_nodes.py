"""Multiplicity of a chosen loop eigenvalue in the degree-four corner characteristic polynomial."""

from collections import defaultdict
from itertools import product
from math import prod
from random import Random

from flint import nmod_mat

from verify_klein_elliptic_completion import rational_example, four_relations
from verify_klein_full_rank_configurations import partner
from verify_klein_two_doublets import modular_points
from verify_quartic_center_order import CentralFiber, integral_quartics

PRIME = 2657


def zero_multiplicity(matrix):
    polynomial = matrix.charpoly()
    return next(i for i in range(polynomial.degree() + 1) if polynomial[i])


def main():
    w, _, reverse = rational_example()
    rules, quartics = integral_quartics(four_relations(w, reverse), PRIME)
    a = modular_points(w, PRIME)
    b = partner(a, w, False, PRIME)
    c = partner(b, reverse, True, PRIME)
    d = partner(c, w.T, False, PRIME)
    rng = Random(8651)
    loop = {word: rng.randrange(PRIME) for word in product(range(4), repeat=4)}
    samples = defaultdict(list)
    for cycle in zip(a, b, c, d):

        def evaluate(section):
            return (
                sum(
                    coeff * prod(row[j] for row, j in zip(cycle, word))
                    for word, coeff in section.items()
                )
                % PRIME
            )

        coords = [evaluate(piece[0]) for piece in quartics]
        scale = pow(coords[0], -1, PRIME)
        samples[tuple(v * scale % PRIME for v in coords)].append(
            evaluate(loop) * scale % PRIME
        )
    assert sorted(map(len, samples.values())) == [1, 1, 1, 1, 4, 4, 4, 4]
    for point, values in samples.items():
        f = CentralFiber(rules, quartics, point[:4], PRIME)
        identity = f.prefix(0, 12, ())
        z = f.fifth_action(quartics[4]) - point[-1] * identity
        u = nmod_mat(64, 64, PRIME)
        for coeff, op in zip(loop.values(), f.loop_operators(loop)):
            u += coeff * op
        assert u * z == z * u and zero_multiplicity(z) == 32
        for value in values:
            delta = u - value * identity
            assert min(zero_multiplicity(delta + c * z) for c in (1, 3, 7)) == 8
        f.prefix.cache_clear()
    zero = CentralFiber(rules, quartics, (0, 0, 0, 0), PRIME)
    assert zero.fifth_action(quartics[4]).minpoly().degree() == 4
    for i in range(4):
        embedding = zero.matrix(
            [zero.evaluate(i, 13, w).entries() for w in zero.basis[i, 9]]
        )
        assert embedding.nrows() == embedding.ncols() and embedding.det()
    generic = CentralFiber(rules, quartics, (1, 2, 3, 4), PRIME)
    polynomial = generic.fifth_action(quartics[4]).minpoly()
    assert (
        polynomial.degree() == 4
        and polynomial.gcd(polynomial.derivative()).degree() == 0
    )
    print(
        "All twenty loop eigenvalues are simple in the reduced degree-four polynomial; central roots have multiplicity two at the sample."
    )


if __name__ == "__main__":
    main()
