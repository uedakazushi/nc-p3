"""Multiplicity-three quartic pencils and polar rank tests over the two split primes."""

from collections import Counter
from itertools import combinations_with_replacement
from math import factorial, prod
from random import Random

from flint import nmod_mat, nmod_poly

from verify_klein_elliptic_completion import rational_example
from verify_klein_full_rank_configurations import partner
from verify_klein_two_doublets import modular_points
from verify_single_configuration_projection import repeated_first_space, square_roots
from verify_unrestricted_doublet_deformations import K, symmetric_example

PAIRS = list(combinations_with_replacement(range(10), 2))
QUARTICS = list(combinations_with_replacement(range(10), 4))


def quartic_pencil(rows, prime):
    index = {monomial: i for i, monomial in enumerate(QUARTICS)}
    entries = []
    for h in rows:
        for rem in PAIRS:
            row = [0] * len(QUARTICS)
            for i, j in PAIRS:
                monomial = tuple(sorted((*rem, i, j)))
                counts = Counter(monomial)
                # Coefficient of s^2 times the fixed quadratic in F(v+s*h).
                multiplicity = (
                    counts[i] * (counts[j] - int(i == j)) // (1 + int(i == j))
                )
                row[index[monomial]] = multiplicity * h[i] * h[j] % prime
            entries.append(row)
    basis, dim = nmod_mat(entries, prime).nullspace()
    assert dim == 2
    return [[int(basis[j, k]) for j in range(len(QUARTICS))] for k in range(dim)]


def polar_matrix(coefficients, point, prime):
    # Symmetric multilinear coefficients, with the common factor 4! omitted.
    symmetric = {
        m: c * prod(factorial(v) for v in Counter(m).values())
        for m, c in zip(QUARTICS, coefficients)
    }
    return nmod_mat(
        [
            [
                sum(
                    point[l] * symmetric[tuple(sorted((i, j, k, l)))] for l in range(10)
                )
                % prime
                for k in range(10)
            ]
            for i, j in PAIRS
        ],
        prime,
    )


def polar_bound(pencil, point, prime):
    polar = [polar_matrix(f, point, prime) for f in pencil]
    interpolate = nmod_mat(
        [[j**i % prime for i in range(8)] for j in range(8)], prime
    ).inv()
    chooser = Random(7301)
    gcd = nmod_poly([], prime)
    infinity = 7
    for _ in range(2):
        # Fixed test projections. Each compressed determinant is a linear
        # combination of true 7-minors by Cauchy--Binet.
        left = nmod_mat(7, 55, [chooser.randrange(prime) for _ in range(7 * 55)], prime)
        right = nmod_mat(
            10, 7, [chooser.randrange(prime) for _ in range(10 * 7)], prime
        )
        a, b = [left * m * right for m in polar]
        values = [int((a + j * b).det()) for j in range(8)]
        coeffs = (interpolate * nmod_mat(8, 1, values, prime)).entries()
        polynomial = nmod_poly([int(v) for v in coeffs], prime)
        assert polynomial
        gcd = gcd.gcd(polynomial)
        infinity = min(infinity, 7 - polynomial.degree())
    assert infinity == 0 and gcd.degree() == 2
    return gcd


def check(points, prime, doubled):
    pairs = list(combinations_with_replacement(K, 2))
    rows = [[q[i] * q[j] % prime for i, j in pairs] for q in points]
    triples = list(combinations_with_replacement(K, 3))
    cubes = [[q[i] * q[j] * q[k] % prime for i, j, k in triples] for q in points]
    assert nmod_mat(cubes, prime).rank() == 20 and nmod_mat(rows, prime).rank() == 10
    pencil = quartic_pencil(rows, prime)
    gcd = polar_bound(pencil, rows[0], prime)
    assert gcd.gcd(gcd.derivative()).degree() == int(doubled)
    square_roots(points, prime, doubled)
    print(
        f"Prime {prime}: 715 coefficients, kernel 2, polar gcd {gcd}, no common root at infinity."
    )


if __name__ == "__main__":
    for w, prime in [(rational_example()[0], 2657), (symmetric_example(), 12619)]:
        if prime == 2657:
            reversed_window = type(w)(
                16, 16, lambda ij, kl: w[4 * (ij // 4) + kl % 4, 4 * (kl // 4) + ij % 4]
            )
            assert repeated_first_space(reversed_window) == [0] * 4
        first = modular_points(w, prime)
        second = partner(first, w, False, prime)
        check(first, prime, prime == 12619)
        check(second, prime, False)
