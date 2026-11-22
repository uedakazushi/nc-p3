"""Exact Gram matrices from two unital ternary spin factors."""
from itertools import combinations_with_replacement

from sympy import Matrix, diag, eye, zeros


def multiply(x, y):
    return Matrix([x.dot(y), x[0]*y[1] + x[1]*y[0], x[0]*y[2] + x[2]*y[0]])


def opposite_rank(g):
    vectors = list(eye(4).columnspace())

    def gram_value(x, y):
        u = multiply(Matrix([x[0], x[2], x[3]]), Matrix([y[1], y[2], y[3]]))
        gy = g*y[2:, :]
        v = multiply(Matrix([x[1], x[2], x[3]]), Matrix([y[0], gy[0], gy[1]]))
        return u.col_join(v)

    values = {(i, j): gram_value(x, y) for i, x in enumerate(vectors)
              for j, y in enumerate(vectors)}
    assert Matrix.hstack(*values.values()).rank() == 6
    h = diag(eye(3), -eye(3)/g.det())
    gram = {(i, j, k, l): (a.T*h*b)[0] for (i, j), a in values.items()
            for (k, l), b in values.items()}
    assert all(value == gram[i, l, k, j] for (i, j, k, l), value in gram.items())
    pairs = list(combinations_with_replacement(range(4), 2))
    w = Matrix([[gram[i, j, k, l] for j, l in pairs] for i, k in pairs])
    assert w[pairs.index((0, 1)), :] == zeros(1, 10)
    assert w[:, pairs.index((0, 1))] == zeros(10, 1)
    leading = [pairs.index(pair) for pair in [(0, 0), (0, 2), (0, 3),
                                              (1, 1), (1, 2), (1, 3)]]
    assert w.extract(leading, leading).rank() == 6
    expected = 6 if g == eye(2)*g[0, 0] else 8
    assert w.rank() == expected
    return expected


def radical_rank(q, alternating):
    # Both axis units coincide with the common radical r. Its relation is
    # the nonsingular binary dot product. The two induced maps are I and
    # q^{-1} q^{-T} and can be aligned or independent.
    vectors = list(eye(4).columnspace())
    h = diag(eye(2), eye(2), Matrix([[0, 1], [1, 0]]))
    omega = Matrix([[0, 1], [-1, 0]])

    def value(x, y):
        u, v = x[2:, :], y[2:, :]
        p = x[0]*v + y[0]*u
        s = x[1]*q*v + y[1]*q.inv().T*u
        height = x[0]*y[0] + x[1]*y[1] + alternating*(u.T*omega*v)[0]
        return p.col_join(s).col_join(Matrix([height, u.dot(v)]))

    values = {(i, j): value(x, y) for i, x in enumerate(vectors)
              for j, y in enumerate(vectors)}
    assert Matrix.hstack(*values.values()).rank() == 6
    gram = {(i, j, k, l): (a.T*h*b)[0] for (i, j), a in values.items()
            for (k, l), b in values.items()}
    assert all(v == gram[i, l, k, j] for (i, j, k, l), v in gram.items())
    pairs = list(combinations_with_replacement(range(4), 2))
    w = Matrix([[gram[i, j, k, l] for j, l in pairs] for i, k in pairs])
    axes = [i for i, (a, b) in enumerate(pairs) if a < 2 and (a, b) != (0, 1)]
    assert w.extract(axes, axes).rank() == 4
    j = q.inv()*q.inv().T
    expected = 6 if not alternating and j == j[0, 0]*eye(2) else 8
    assert w.rank() == expected
    return expected


if __name__ == "__main__":
    tests = [eye(2), 3*eye(2), diag(2, 3), Matrix([[1, 2], [0, 1]]),
             Matrix([[1, 3], [2, 4]])]
    print("Transverse factor ranks over Q:", [opposite_rank(g) for g in tests])
    print("Radical-line ranks over Q:", [radical_rank(q, a) for q, a in
          [(eye(2), 0), (eye(2), 1), (diag(1, 2), 0), (diag(1, 2), 1)]])
