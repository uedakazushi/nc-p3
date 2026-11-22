"""Word reversal, half-turn, character translation and an oriented cubic invariant."""

from itertools import product
import sympy as sp
from verify_klein_elliptic_completion import gram_twist
from verify_fixed_exchange_isogeny import exchange
from verify_rank_fifteen_adjoint_invariants import (
    G,
    N,
    P,
    basis,
    mat,
    modp,
    opposite_matrix,
    trace,
    traceless_map,
)


def prove_tensor_identities():
    a, b, e, f, z, lam, u = sp.symbols("a b e f z lambda u")
    r, d, q = sp.symbols("r d q", nonzero=True)
    ii = sp.I

    def tensor(a, b, e, f):
        left = sp.Matrix([[a, b], [z * b, a]])
        right = sp.Matrix([[e, f], [z * f, e]])
        return gram_twist(left, right, z, lam)[0]

    original = tensor(a, b, e, f)
    swapped = tensor(z * b, a, z * f, e)
    reduction = exchange(1, z, lam) * sp.Matrix([a, b])
    change = {e: u * reduction[0], f: u * reduction[1]}

    def verify(target, indices, diagonals, *, use_exchange=False):
        for word in product(range(4), repeat=4):
            i, j, k, t = word
            I, J, K, T = indices(*word)
            if i ^ j ^ k ^ t:
                assert (
                    target[4 * i + j, 4 * k + t] == original[4 * I + J, 4 * K + T] == 0
                )
                continue
            scale = sp.sympify(sp.prod(diagonals[p][letter] for p, letter in enumerate(word)))
            powers = scale.as_powers_dict()
            for root, square in ((r, z), (d, lam), (q, u * (z * b * b - a * a))):
                exponent = powers[root]
                assert exponent % 2 == 0
                scale = scale / root**exponent * square ** (exponent // 2)
            diff = target[4 * i + j, 4 * k + t] - scale * original[4 * I + J, 4 * K + T]
            if use_exchange:
                diff = diff.subs(change)
            assert sp.cancel(diff) == 0, word

    verify(
        swapped,
        lambda i, j, k, t: (t, k, j ^ 1, i ^ 1),
        (
            (1 / d, 1 / (r * d), 1 / r, 1),
            (1 / d, r / d, 1 / r, 1),
            (d, d / r, 1 / r, 1),
            (z * d, r * d, r, 1),
        ),
    )
    perm = (0, 1, 3, 2)
    verify(
        original,
        lambda i, j, k, t: (perm[k], perm[t], perm[i], perm[j]),
        (
            (r / d, r / d, 1, z),
            (1 / (r * d), r / d, 1, 1),
            (d / r, d / r, 1 / z, 1),
            (r * d, d / r, 1, 1),
        ),
    )
    verify(
        original,
        lambda i, j, k, t: (i ^ 2, j ^ 3, k ^ 2, t ^ 3),
        (
            (
                1 / (lam * q),
                1 / (lam * q * r),
                ii / (lam * q * q * r),
                ii / (lam * q * q),
            ),
            (-ii * q, -ii * q * r, 1 / r, 1),
            (-ii * lam * q, -ii * lam * q / r, 1 / r, 1),
            (lam * q * z, lam * q * r, ii * r, ii),
        ),
        use_exchange=True,
    )


def oriented_value(opposite):
    T = modp(traceless_map(opposite))
    gram = modp(G)
    S = gram.inv() * T.transpose() * gram * T
    base = [modp(x) for x in basis]
    dual_image = S * gram.inv()
    images = []
    for i in range(15):
        column = mat(4, 4, [0] * 16, P)
        for j in range(15):
            column += base[j] * dual_image[j, i]
        images.append(column)
    total = 0
    for i, j in product(range(15), repeat=2):
        first = base[i] * base[j] + base[j] * base[i]
        second = images[i] * images[j] + images[j] * images[i]
        for k in range(15):
            total += trace(first * base[k]) * trace(second * images[k])
    return int(total / trace(S) ** 3)


def one_shift(M):
    # New W(i,j,k,l) = old W(l,i,j,k).
    return mat(
        16, 16, [M[4 * t + j, 4 * i + k] for i, k, j, t in product(range(4), repeat=4)]
    )


if __name__ == "__main__":
    prove_tensor_identities()
    for color in range(2):
        M = opposite_matrix(4, 3 * pow(4, -1, N) % N, color)
        pair = oriented_value(M), oriented_value(one_shift(M))
        assert pair == ((10, 44), (44, 10))[color], pair
        print(
            "Presentation", color, "anchored and shifted cubic invariants mod 101", pair
        )
