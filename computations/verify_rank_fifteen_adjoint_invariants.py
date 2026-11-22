"""Two independent indexed-isomorphism invariants on each rank-fifteen surface."""

from itertools import product
import sympy as sp
from flint import nmod_mat, nmod
from verify_klein_elliptic_completion import gram_twist

P = 101
N = P * P


def mat(rows, cols, vals, m=N):
    return nmod_mat(rows, cols, [int(a) % m for a in vals], m)


def matrix_unit(i, j):
    return mat(4, 4, [int((i, j) == (a, b)) for a, b in product(range(4), repeat=2)])


basis = [matrix_unit(i, j) for i, j in product(range(4), repeat=2) if i != j] + [
    matrix_unit(i, i) - matrix_unit(3, 3) for i in range(3)
]


def trace(x):
    return sum(x[j, j] for j in range(x.nrows()))


G = mat(15, 15, [trace(a * b) for a in basis for b in basis])


def modp(A):
    return mat(A.nrows(), A.ncols(), A.entries(), P)


def inverse_mod_prime_squared(A):
    """Lift the inverse from the prime field by one Newton step."""
    d = A.nrows()
    B = mat(d, d, modp(A).inv().entries())
    eye = mat(d, d, [int(i == j) for i, j in product(range(d), repeat=2)])
    ret = B * (eye * 2 - A * B)
    assert A * ret == eye
    return ret


def kernel_form(A):
    """Recover the null form from a rank-15 minor with unit determinant."""
    d = 16
    D = modp(A)
    R, r = D.nullspace()
    L, l = D.transpose().nullspace()
    assert l == r == 1
    row = next(i for i in range(d) if L[i, 0])
    col = next(i for i in range(d) if R[i, 0])
    cols = [i for i in range(d) if i != col]
    rows = [i for i in range(d) if i != row]
    U = mat(15, 15, [A[i, j] for i in rows for j in cols])
    q = inverse_mod_prime_squared(U) * mat(15, 1, [-A[i, col] for i in rows])
    v = [q[cols.index(i), 0] if i != col else 1 for i in range(d)]
    assert not any((A * mat(16, 1, v)).entries())
    t = mat(4, 4, v)
    assert modp(t).det()
    return t


def lift_root(val, root):
    root = nmod(root, N)
    assert int(root * root - val) % P == 0
    root = (root + val / root) / 2
    assert root * root == val
    return root


def opposite_matrix(x, c, color):
    x = nmod(x, N)
    c = nmod(c, N)
    L = lift_root(1 + c * c, 5 * pow(4, -1, P) % P)
    s = c + 1 / c
    k = (x + 1 / x) / 2
    H = (x - 1 / x) / 2
    n = H / L
    A = lift_root((k + s) ** 2 - 1, 10)
    B = lift_root((k - s) ** 2 - 1, 6)
    r = k + (H * H + H * A) / s
    w = k + (-H * H + H * B) / s
    z = r * r
    l = w / r
    assert int(z * (z * z - 1) * (w * w - 1) * (w * w - z) * (w * w * z - 1)) % P
    a, b, e, f = 1 + n, (k + c * n) / r, 1 - n, (k - c * n) / r
    if color:
        a, b, e, f = z * b, a, z * f, e
    z, l, a, b, e, f = map(lambda v: sp.Integer(int(v)), (z, l, a, b, e, f))
    S = sp.Matrix([[a, b], [z * b, a]])
    T = sp.Matrix([[e, f], [z * f, e]])
    _, opp, _ = gram_twist(S, T, z, l)
    return mat(16, 16, [int(q.p) * pow(int(q.q), -1, N) for q in opp])


def invariants(x, c, color):
    M = opposite_matrix(x, c, color)
    values = opposite_invariants(M)
    changes = []
    for i in range(4):
        change = sp.eye(4)
        change[i, (i + 1) % 4] = i + 1
        changes.append(change)
    rows = mat(16, 16, sp.kronecker_product(changes[0], changes[2]))
    columns = mat(16, 16, sp.kronecker_product(changes[1], changes[3]))
    assert opposite_invariants(7 * rows * M * columns.transpose()) == values
    return values


def traceless_map(M):
    right = kernel_form(M)
    left = kernel_form(M.transpose())
    cols = []
    # The input identity maps to the right null form. The output uses the left form as trace.
    for X in basis:
        raw = X.transpose() * right
        vec = M * mat(16, 1, raw.entries())
        Y = mat(4, 4, vec.entries()) * left.transpose()
        assert trace(Y) == 0
        cols.append(
            [Y[i, j] for i, j in product(range(4), repeat=2) if i != j]
            + [Y[i, i] for i in range(3)]
        )
    T = mat(15, 15, [v for col in cols for v in col]).transpose()
    assert modp(T).det()
    return T


def opposite_invariants(M):
    T = traceless_map(M)
    S = inverse_mod_prime_squared(G) * T.transpose() * G * T
    assert int(trace(S)) % P
    return tuple(int(trace(S**m) / trace(S) ** m) for m in (2, 3))


if __name__ == "__main__":
    for co in range(2):
        x = 4
        c = 3 * pow(4, -1, N) % N
        v = invariants(x, c, co)
        d = []
        for u in (invariants(x + P, c, co), invariants(x, c + P, co)):
            assert all((b - a) % P == 0 for a, b in zip(v, u))
            d.append([(b - a) // P % P for a, b in zip(v, u)])
        determinant = (d[0][0] * d[1][1] - d[1][0] * d[0][1]) % P
        assert v == (1418, 4898) and d == [[22, 88], [23, 66]] and determinant == 34
        print(
            "Presentation",
            co,
            "Jacobian columns",
            d,
            "determinant modulo 101",
            determinant,
        )
