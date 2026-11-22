"""Exact separated block ranks, including an everywhere singular square pencil."""
from itertools import combinations_with_replacement

from sympy import Matrix, diag, eye, zeros


def kernel_dimension(a, b):
    s = a.rows
    z = zeros(s, 4)
    system = Matrix.vstack(Matrix.hstack(a,z,b,z), Matrix.hstack(z,a,z,b),
                           Matrix.hstack(a,b,z,z), Matrix.hstack(z,z,a,b))
    return system.cols-system.rank()


if __name__ == "__main__":
    nilpotent = Matrix([[0,1,0,0], [0,0,1,0], [0,0,0,1], [0,0,0,0]])
    cases = [(eye(4), nilpotent), (nilpotent, eye(4)),
             (diag(1,1,1,0), diag(0,1,1,0)), (zeros(4), zeros(4))]
    values = [kernel_dimension(a,b) for a,b in cases]
    assert values == [4, 4, 7, 16]
    print("Separated block dimensions over Q:", values)
    for r in (3, 4):
        evaluation = zeros(6, 4)
        evaluation[:r, :r] = eye(r)
        columns = []
        for i, j in combinations_with_replacement(range(4), 2):
            q = zeros(4)
            q[i, j] = q[j, i] = 1
            columns.append((evaluation*q).reshape(24, 1))
        last = Matrix.hstack(*columns)
        assert last.rank() == (9 if r == 3 else 10)
        if r == 3:
            assert last.nullspace() == [eye(10)[:, -1]]
    print("Last Eagon-Northcott ranks 9 and 10; kernel is the Veronese square.")
