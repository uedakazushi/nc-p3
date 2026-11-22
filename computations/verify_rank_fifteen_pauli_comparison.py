"""Finite indexed-isomorphism comparison in a paired four-dimensional Pauli basis."""

from itertools import combinations, product
import numpy as np
from flint import nmod
from verify_rank_fifteen_endpoint_orbits import trace_blocks
from verify_rank_fifteen_word_symmetries import prove_tensor_identities
from verify_rank_fourteen_geometric_reduction import verify_central_profile
from rank_fifteen_octic import finite_windows


def tensor_operator(frame, p):
    """The opposite map in End(V0) tensor End(V1), viewed on V0 tensor V1."""
    T, _ = trace_blocks(frame)
    M = np.zeros((4, 4, 4, 4), dtype=np.int64)
    for h, i, j in product(range(4), repeat=3):
        M[i, j ^ h, i ^ h, j] = int(T[h][i][j])
    return M.reshape(16, 16)


def swap(M):
    return M.reshape(4, 4, 4, 4).transpose(0, 3, 2, 1).reshape(16, 16)


def closure(M, p):
    """Span generated from I,H by multiplication, also reversing the second factor."""
    I = np.eye(16, dtype=np.int64)
    base = []
    piv = {}

    def add(X):
        line = X.flatten()
        for j, row in piv.items():
            if line[j]:
                line = (line - line[j] * row) % p
        hits = np.flatnonzero(line)
        if len(hits):
            j = hits[0]
            piv[j] = line * pow(int(line[j]), -1, p) % p
            base.append(X)

    add(I)
    add(M)
    i = 0
    while i < len(base):
        for j in range(i + 1):
            for a, b in [(base[j], base[i]), (base[i], base[j])]:
                add(a @ b % p)
                add(swap(swap(a) @ swap(b) % p))
        i += 1
    return len(base)


def weighted(vals, shift):
    return np.array(
        [
            [int(v) if j == i ^ shift else 0 for j in range(4)]
            for i, v in enumerate(vals)
        ],
        dtype=np.int64,
    )


def form(x, y, side):
    a, b = x % 4, x // 4
    c, d = y % 4, y // 4
    if side:
        b, d = (0, 3, 2, 1)[b], (0, 3, 2, 1)[d]
    return ((a & d).bit_count() + (b & c).bit_count()) % 2


def label_permutations():
    """Isometries of both alternating commutator forms on the four-bit labels."""
    basis = (1, 2, 4, 8)
    for ims in product(range(1, 16), repeat=4):
        if all(
            form(ims[i], ims[j], side) == form(basis[i], basis[j], side)
            for i, j in combinations(range(4), 2)
            for side in (0, 1)
        ):
            F = [0]
            for k in ims:
                F += [v ^ k for v in F]
            yield F


def table(frame, r, d, p):
    """Recover all fifteen nonidentity terms, then multiply their normalized representatives."""
    z = frame[0]
    assert r * r == z and d * d == frame[1]
    A = [weighted((r, r, 1, z), 1), weighted((1 / r, r, 1, 1), 1)]
    B = [weighted((1, 1 / r, d, r * d), 2), weighted((1, r, d, r * d), 3)]
    I = np.eye(4, dtype=np.int64)
    big = np.eye(16, dtype=np.int64)
    for x, y in zip(A, B):
        assert np.array_equal(x @ y % p, y @ x % p)
    H = tensor_operator(frame, p)
    basis = []
    for l in range(16):
        a, b = l % 4, l // 4
        pair = []
        for j in range(2):
            shift = (0, 1, 3, 2)[a] if j else a
            phase = (0, 3, 2, 1)[b] if j else b
            Z = np.diag([(-1) ** ((phase & k).bit_count()) for k in range(4)])
            pair.append((A[j] if shift & 1 else I) @ (B[j] if shift & 2 else I) @ Z % p)
        G = np.kron(*pair) % p
        scale = int((G @ G)[0, 0]) % p
        assert np.array_equal((G @ G) % p, scale * big % p)
        c = int((G @ H).trace()) * pow(16 * scale, -1, p) % p
        assert bool(c) == bool(l)
        basis.append(G * c % p if l else big)
    assert np.array_equal(sum(basis[1:]) % p, H)
    D = np.zeros((16, 16), dtype=np.int64)
    for x, y in product(range(16), repeat=2):
        K = basis[x ^ y]
        i, j = np.argwhere(K)[0]
        val = basis[x] @ basis[y] % p
        D[x, y] = int(val[i, j]) * pow(int(K[i, j]), -1, p) % p
        assert np.array_equal(val, D[x, y] * K % p)
    return D


def match(A, B, perms, p):
    """A single scalar is forced by the distinct labels 1 and 2."""
    ans = []
    for F in perms:
        C = B[np.ix_(F, F)]
        a = int(C[1, 2]) * pow(int(A[1, 2]), -1, p) % p
        E = np.full((16, 16), a, dtype=np.int64)
        np.fill_diagonal(E, a * a % p)
        E[0, :] = E[:, 0] = 1
        if np.array_equal(A * E % p, C):
            ans.append(F)
    return ans


def verify():
    prove_tensor_identities()
    prime = 101
    r, d = nmod(2, prime), nmod(3, prime)
    z, lam = r * r, d * d
    k = (lam + 1) * (lam * z + 1) / (4 * lam * z)
    m = -(lam - 1) * (lam * z - 1) / (4 * lam * z)
    labels = list(label_permutations())
    assert len(labels) == 48
    roots = [
        nmod(n, prime)
        for n in range(prime)
        if not (nmod(n, prime) ** 4 + (1 - z * k * k) * n * n + z * m * m)
    ]
    assert list(map(int, roots)) == [46, 47, 54, 55] and z * k * k != 1
    for n in roots:
        frame = (z, lam, 1 + n, k + m / n, 1 - n, k - m / n)
        assert all(frame) and closure(tensor_operator(frame, prime), prime) == 16
        unit = nmod(1, prime)
        assert verify_central_profile(
            finite_windows(frame, unit, lambda value: value), unit
        )
    n = roots[0]
    a, b, e, f = 1 + n, k + m / n, 1 - n, k - m / n
    frame = (z, lam, a, b, e, f)
    original = table(frame, r, d, prime)
    reciprocal = table((1 / z, 1 / lam, a, z * b, e, z * f), 1 / r, 1 / d, prime)
    other = table((z, 1 / (lam * z), a, b, e, f), r, 1 / (r * d), prime)
    assert len(match(original, original, labels, prime)) == 1
    assert len(match(original, reciprocal, labels, prime)) == 1
    assert len(match(original, table(frame, -r, -d, prime), labels, prime)) == 1
    assert not match(original, other, labels, prime)
    print(
        "Pauli algebra dimension 16; joint isometry group 48; one generic stabilizer label action. All four interior residue roots have full central cuts."
    )


if __name__ == "__main__":
    verify()
