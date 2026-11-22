"""Recover regular orbit coordinates from projective tensor contractions at every point."""

from collections import Counter
from itertools import product

import numpy as np
import sympy as sp
from flint import nmod

from verify_rank_fifteen_coordinate_collisions import oriented_cubic
from verify_rank_fifteen_endpoint_orbits import multiply, trace, trace_blocks
from verify_rank_fifteen_parameter_orbits import dot, pattern, product_sign, quadrics
from verify_rank_fifteen_pauli_comparison import form, tensor_operator
from verify_rank_fifteen_rational_quotient import relative_forms
from verify_rank_fifteen_root_exchange import coefficients


def cycle_sign(labels, side):
    prefix = 0
    parity = 0
    for label in labels:
        phase = prefix // 4
        if side:
            phase = (0, 3, 2, 1)[phase]
        parity += dot(phase, label % 4)
        prefix ^= label
    assert prefix == 0
    return parity


def odd_fifth_terms():
    """The right 3-cycle and 2-cycle force labels a,b,a+b and c,c."""
    terms = Counter()
    for a, b, c in product(range(1, 16), repeat=3):
        if a == b:
            continue
        labels = (a, c, b, a ^ b, c)
        if sum(dot(j % 4, j // 4) for j in labels) % 2 == 0:
            continue
        sign = cycle_sign([labels[j] for j in (0, 4, 2, 3, 1)], 0)
        sign += cycle_sign([labels[j] for j in (0, 2, 3)], 1)
        sign += cycle_sign([labels[j] for j in (1, 4)], 1)
        terms[tuple(sorted(labels))] += (-1) ** sign
    return terms


def verify_symbolic():
    variables = sp.symbols("a b L M A B C D")
    row = pattern(*variables)
    p2, p3, d3, h, E = relative_forms(*variables)
    ideal = sp.groebner(
        quadrics(*variables), *variables, domain=sp.QQ_I, order="grevlex"
    )

    def zero(expr):
        assert ideal.reduce(sp.expand(expr))[1] == 0

    mu = [(-1) ** product_sign(j, j) * v * v for j, v in enumerate(row)]
    zero(sum(mu) - p2)
    zero(4 * sum(z * z for z in mu) - p2 * p2 + 3 * h * h)
    cube = 0
    difference = 0
    for a, b in product(range(1, 16), repeat=2):
        c = a ^ b
        if c == 0:
            continue
        sign = cycle_sign((a, b, c), 0) + cycle_sign((a, b, c), 1)
        cube += (-1) ** sign * row[a] * row[b] * row[c]
        difference += (form(a, b, 0) - form(a, b, 1)) * mu[a] * mu[b] * mu[c]
    zero(cube - 6 * p3)
    zero(4 * difference - 3 * E)
    odd = sum(k * sp.prod(row[j] for j in js) for js, k in odd_fifth_terms().items())
    zero(odd + h * d3)


def verify_frames():
    prime = 101
    unit = nmod(1, prime)
    for z, lam in ((4, 9), (4, 87), (21, 45), (21, 87)):
        z, lam, n = z * unit, lam * unit, 46 * unit
        k, m = coefficients(z, lam)
        U, V = z * k * k, z * m * m
        frame = (z, lam, 1 + n, k + m / n, 1 - n, k - m / n)
        tensor = tensor_operator(frame, prime)

        def trace_power(exponent):
            matrix = tensor
            for _ in range(1, exponent):
                matrix = matrix @ tensor % prime
            return unit * int(matrix.trace())

        p2, p3 = trace_power(2) / 16, trace_power(3) / 96
        blocks, inputs = trace_blocks(frame)
        outputs = [
            multiply(A, [[A[j ^ h][i ^ h] for j in range(4)] for i in range(4)])
            for h, A in enumerate(blocks)
        ]
        h2 = (p2 * p2 - sum(trace(multiply(S, S)) for S in inputs) / 64) / 3
        assert (p2, p3, h2) == (U + 3, U - 1, V)
        q = (oriented_cubic(inputs) - oriented_cubic(outputs)) / (3072 * p3 * p3)
        assert q == (2 * n * n + 1 - U) * (z - 1) * (z * lam * lam - 1) / (
            4 * z * lam * (U - 1) ** 2
        )
        P = (4, 0, 3, 1, 2)
        Q = (2, 4, 3, 0, 1)
        labels = "abcdefghij"
        contraction = (
            ",".join(
                labels[j] + labels[5 + j] + labels[P[j]] + labels[5 + Q[j]]
                for j in range(5)
            )
            + "->"
        )
        assert 4**10 * (prime - 1) ** 5 < 2**63

        def fifth(A):
            return unit * int(np.einsum(contraction, *([A] * 5), optimize=True))

        A = tensor.reshape(4, 4, 4, 4)
        p = -h2 * (fifth(A) - fifth(A.transpose(2, 1, 0, 3))) / (128 * p3**3)
        assert p == V * z * k * m * (U - 9) / (U - 1) ** 3


if __name__ == "__main__":
    verify_symbolic()
    verify_frames()
    print(
        "All three projective quotient coordinates are intrinsic on every ordinary rank-fifteen window"
    )
