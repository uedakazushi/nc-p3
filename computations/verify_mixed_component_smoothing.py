"""Schur equations at the split mixed window and its transverse polynomial curve."""
import sympy as sp
import numpy as np
from verify_semisimple_period_obstruction import FRAME, G, P, Q, configurations, reduced


def window(g, f, q, p):
    tensor = np.zeros((4, 4, 4, 4), dtype=np.int64)
    for i in range(3):
        tensor[0, i+1, 0, i+1] = tensor[i+1, 0, i+1, 0] = 1
        tensor[i+1, 1:, 0, 1:] = g[i]
        tensor[0, 1:, i+1, 1:] = f[i]
        tensor[1:, i+1, 1:, 0] = q[i]
        tensor[1:, 0, 1:, i+1] = p[i]
    tensor[1:, 1:, 1:, 1:] = np.einsum('abe,ced->abcd', g, f) + np.einsum('bae,dec->abcd', q, p)
    return tensor


def along(t):
    g, f, q, p = G+t*Q, -G+t*Q, Q+t*G, Q-t*G
    for a, c, b, d in np.ndindex(3, 3, 3, 3):
        assert (g[a]@f[c]-f[c]@g[a])[b, d] + (q[b]@p[d]-p[d]@q[b])[a, c] == 0
    return window(g, f, q, p)


def kernel(matrix):
    array, pivots = reduced(matrix)
    free = [i for i in range(matrix.shape[1]) if i not in pivots]
    basis = np.zeros((matrix.shape[1], len(free)), dtype=np.int64)
    basis[free] = np.eye(len(free), dtype=np.int64)
    basis[pivots] = -array[:len(pivots), free]
    return basis.T % P


def rank(matrix):
    return len(reduced(matrix)[1])


def main():
    w = along(0)
    gradients = []
    for axes in [(0, 1, 2, 3), (2, 1, 0, 3)]:
        matrix = w.transpose(axes).reshape(16, 16)
        left, right = kernel(matrix.T), kernel(matrix)
        assert len(left) == len(right) == 10
        gradients.extend(np.outer(x, y).reshape(4, 4, 4, 4).transpose(axes).flatten()
                         for x in left for y in right)
    assert rank(gradients) == 180
    odd, a, c = configurations()
    rows = []
    for x, b, y in zip(a, odd, c):
        left = np.kron(FRAME, b[:, None]).T @ matrix
        right = matrix @ np.kron(x[:, None], FRAME)
        assert rank(left) == rank(right) == 3
        assert not np.any(y @ left % P) and not np.any(right @ b % P)
        rows.append(np.kron(np.kron(y, b), np.kron(x, b)))
    assert rank(rows) == 20

    def opposite(w):
        return w.transpose(0, 2, 1, 3).reshape(16, 16)

    base, plus, minus = [opposite(along(t)) for t in (0, 1, -1)]
    assert rank(kernel(base.T) @ (plus-minus) @ kernel(base).T) == 6
    t = sp.symbols('t')
    mat = sp.Matrix(base) + t*sp.Matrix(plus-minus)/2 + t*t*sp.Matrix(plus+minus-2*base)/2
    for color in range(4):
        indices = [4*i + (i ^ color) for i in range(4)]
        other = [i for i in range(16) if i not in indices]
        assert mat.extract(indices, other).is_zero_matrix
        det = mat.extract(indices, indices).det()
        assert sp.expand(det - (16*t*t if color else -12*(1+t*t)**2)) == 0
    # At nonzero c the four diagonal rescalings put this curve in Walsh form.
    c = sp.symbols('c')
    relation = c*c + 2*(1+t*t)
    phases = [[sp.I, 1, sp.I, 1], [sp.I, sp.I, 1, 1],
              [-sp.I, 1, -sp.I, 1], [-sp.I, -sp.I, 1, 1]]
    signed_squares = sp.zeros(3)
    for b, e in np.ndindex(4, 4):
        d = b ^ e
        entries = []
        for i in range(4):
            labels = (i, i^b, i^d, i^e)
            scale = sp.prod(phases[h][j]*(c if j == 0 else 1) for h, j in enumerate(labels))
            x = scale*mat[4*i+(i^d), 4*(i^b)+(i^e)] * (-1)**((i >> 1)*(d & 1))
            entries.append(sp.rem(x, relation, c))
        assert all(sp.expand(x-entries[0]) == 0 for x in entries)
        if not (b and e):
            assert entries[0] == 0
        else:
            square = (-1)**d.bit_count() * entries[0]**2
            signed_squares[b-1, e-1] = sp.cancel(sp.rem(square, relation, c) / (-2*(1+t*t)))
    row = [-2*(1+t*t), (1-t)**2, (1+t)**2]
    circulant = sp.Matrix(3, 3, lambda i, j: row[(j-i) % 3])
    assert (signed_squares-circulant).applyfunc(sp.cancel).is_zero_matrix
    assert sp.factor(circulant[:2, :2].det()) == (t*t+3)*(3*t*t+1)
    print("PASS: transverse rank 6, Schur rank 180, 20 conormals, Walsh normalization and AS exceptions")


if __name__ == '__main__':
    main()
