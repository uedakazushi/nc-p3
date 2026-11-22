"""The simultaneous character flip, its two fixed lines, and two frames at a regular point."""

from itertools import product

import sympy as sp

from verify_klein_elliptic_completion import coeffs, four_relations, gram_twist
from verify_klein_period_four_center import artinian_quotient, commuting_space, paths


def flip(matrix, z):
    a, b, c, d = matrix
    return sp.Matrix([[d, c / z], [z * b, a]])


def check_action():
    q, l = sp.symbols("q l", nonzero=True)
    z = q * q
    a, b, c, d, e, f, g, h = sp.symbols("a b c d e f g h")
    A, B = sp.Matrix([[a, b], [c, d]]), sp.Matrix([[e, f], [g, h]])
    # The thick blocks are determined by their upper left and two cross blocks.
    J = sp.Matrix([[0, 1], [1, 0]])
    low = sp.diag(1, z)
    high = sp.diag(1, z)
    X = sp.Matrix([[1, 1 / l], [1 / l, z]])
    Y = sp.Matrix([[1, l], [l * z, 1]])
    assert (low * J * X * J * high / z - X).applyfunc(sp.cancel) == sp.zeros(2)
    assert (low * J * Y * J * low.inv() - Y).applyfunc(sp.cancel) == sp.zeros(2)
    for M in (A, B):
        assert low * J * M * J * low.inv() == flip(M, z)
    alpha = (q, q, 1, z)
    beta = (1 / z, 1, 1 / q, 1 / q)
    gamma = alpha
    delta = (1, 1 / z, 1 / q, 1 / q)

    def thin(i, j):
        return z if (i, j) == (3, 1) else 1

    for i, j, k, n in product(range(4), repeat=4):
        if i ^ j ^ k ^ n or (i ^ j) < 2:
            continue
        assert (
            sp.cancel(
                alpha[i]
                * beta[j]
                * gamma[k]
                * delta[n]
                * thin(i ^ 1, j ^ 1)
                * thin(k ^ 1, n)
                - thin(i, j) * thin(k, n ^ 1)
            )
            == 0
        )
    for sign in (1, -1):
        co, _, _ = coeffs(a, b, sign * z * b, sign * a, z, l)
        assert (
            sp.cancel(co[3] - sign * z * co[1]) == sp.cancel(co[7] - sign * co[5]) == 0
        )


def regular_example():
    A = sp.Matrix([[0, 1], [3, 0]])
    B = sp.Matrix([[2, sp.Rational(8, 13)], [sp.Rational(24, 13), 2]])
    w, opp, rev = gram_twist(A, B, sp.Integer(3), sp.Integer(2))
    assert (w.rank(), opp.rank(), rev.rank()) == (6, 16, 6)
    rel = four_relations(w, rev)
    basis, image = paths(rel, 5)
    assert [commuting_space(basis, image, h, 0).cols for h in range(4)] == [3, 2, 0, 0]
    artinian_quotient(rel, basis, image)
    prime = 1009
    q = int(sp.sqrt_mod(3, prime))
    im = int(sp.sqrt_mod(-1, prime))
    invq = pow(q, -1, prime)
    weights = (
        (q, q, 1, 3),
        (pow(3, -1, prime), 1, invq, invq),
        (q, q, 1, 3),
        (1, pow(3, -1, prime), invq, invq),
    )
    changes = []
    for j, weight in enumerate(weights):
        root = im * (q if j % 2 == 0 else invq) % prime
        basis = sp.diag(
            *[sp.Matrix([[weight[i], weight[i]], [root, -root]]) for i in (0, 2)]
        )
        changes.append(basis.inv_mod(prime))
    indices = list(product(range(4), repeat=4))
    entries = {key: w[4 * key[0] + key[1], 4 * key[2] + key[3]] for key in indices}
    entries = {
        key: int(v.p) * pow(int(v.q), -1, prime) % prime for key, v in entries.items()
    }
    for j, T in enumerate(changes):
        entries = {
            key: sum(
                T[key[j], i] * entries[key[:j] + (i,) + key[j + 1 :]] for i in range(4)
            )
            % prime
            for key in indices
        }
    assert all(v == 0 for (i, j, k, n), v in entries.items() if i ^ j ^ k ^ n)
    for reverse, h in product(range(2), range(4)):
        block = sp.Matrix(
            4,
            4,
            lambda i, k: entries[
                (k, i ^ h, i, k ^ h) if reverse else (i, i ^ h, k, k ^ h)
            ],
        )
        if h < 2:
            assert all(
                block[i : i + 2, j : j + 2].det() % prime
                for i, j in product((0, 2), repeat=2)
            )
        else:
            assert all(block)
    print(
        "The full-rank AS window has a second genuinely distinct marking inside the invertible chart."
    )


if __name__ == "__main__":
    check_action()
    regular_example()
