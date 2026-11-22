"""Exact orbit-normal tangents and second-order tests for two four-windows."""

from itertools import combinations, combinations_with_replacement, product

import sympy as sp

from verify_klein_elliptic_completion import rational_example
from verify_klein_two_doublets import cross_blocks, gram, other_completion

K = range(4)


def rank(matrix):
    return matrix.to_DM().rank()


def independent(matrix):
    return matrix[:, list(matrix.to_DM().rref()[1])]


def kernel(matrix):
    rows = matrix.to_DM().nullspace().to_Matrix()
    return [rows[i, :].T / next(x for x in rows[i, :] if x) for i in range(rows.rows)]


def exchange(matrix):
    """Interchange the first and third letters before taking the adjacent cut."""
    return sp.Matrix(
        16, 16, lambda ij, kl: matrix[4 * (kl // 4) + ij % 4, 4 * (ij // 4) + kl % 4]
    )


def opposite_cut(matrix):
    return sp.Matrix(
        16, 16, lambda ik, jl: matrix[4 * (ik // 4) + jl // 4, 4 * (ik % 4) + jl % 4]
    )


def normal_data(matrix):
    left = sp.Matrix.hstack(*kernel(matrix.T))
    right = sp.Matrix.hstack(*kernel(matrix))
    cols = list(matrix.to_DM().rref()[1])
    rows = list(matrix.T.to_DM().rref()[1])
    pivot_inverse = matrix[rows, cols].inv()
    inverse = sp.zeros(16)
    for i, j in product(range(len(cols)), repeat=2):
        inverse[cols[i], rows[j]] = pivot_inverse[i, j]
    assert matrix * inverse * matrix == matrix

    def project(tensor):
        normal = left.T * tensor * right
        return normal.reshape(len(normal), 1)

    return project, inverse


def orbit(w, h, sign):
    columns = []
    for i in K:
        e = sp.zeros(4)
        e[i, i ^ h] = 1
        even = sp.kronecker_product(e, sp.eye(4))
        odd = sp.kronecker_product(sp.eye(4), e)
        if sign == 0:
            columns.extend([even * w, w * even.T, odd * w, w * odd.T])
        else:
            columns.append(even * w + sign * w * even.T)
            if sign == 1:
                columns.extend([odd * w, w * odd.T])
    return independent(sp.Matrix.hstack(*[t.reshape(256, 1) for t in columns]))


def tensor_basis(h, sign):
    if sign == 0:
        pairs = product(K, repeat=2)
    elif sign == 1:
        pairs = combinations_with_replacement(K, 2)
    else:
        pairs = combinations(K, 2)
    for i, k in pairs:
        for j in K:
            l = h ^ i ^ j ^ k
            term = sp.zeros(16)
            term[4 * i + j, 4 * k + l] = 1
            if sign:
                term[4 * k + j, 4 * i + l] = sign
            yield term


def tangent(w, h, sign, project):
    basis = list(tensor_basis(h, sign))
    derivative = sp.Matrix.hstack(*map(project, basis))
    change = orbit(w, h, sign)
    assert all(
        project(change[:, i].reshape(16, 16)).is_zero_matrix for i in range(change.cols)
    )
    orbit_dim = change.cols
    lifts = []
    for coeff in kernel(derivative):
        tensor = sum((c * b for c, b in zip(coeff, basis) if c), sp.zeros(16))
        enlarged = change.row_join(tensor.reshape(256, 1))
        if rank(enlarged) > change.cols:
            change = enlarged
            lifts.append(tensor)
    print(
        f"sign {sign}, character {h}: orbit {orbit_dim}, quotient tangent {len(lifts)}",
        flush=True,
    )
    return derivative, lifts, orbit_dim


def symmetric_example():
    u = sp.Matrix([[1, 2], [3, 5]])
    v = sp.Matrix([[2, 3], [5, 7]])
    a = sp.Matrix([[2, 1], [5, 3]])
    return gram(u, v, a, other_completion(a, *cross_blocks(u, v)))


def rank_ten():
    w = symmetric_example()
    assert w == exchange(w)
    project, inv = normal_data(w)
    plus = [tangent(w, h, 1, project) for h in K]
    minus = [tangent(w, h, -1, project) for h in K]
    assert sum(rank(mat) for mat, _, _ in plus) == 100
    assert [len(t) for _, t, _ in plus] == [4, 2, 4, 4]
    assert [len(t) for _, t, _ in minus] == [1, 0, 0, 0]
    assert sum(dim for _, _, dim in plus + minus) == 61
    (anti,) = minus[0][1]
    for h, (_, lifts, _) in enumerate(plus):
        mixed = sp.Matrix.hstack(
            *[project(t * inv * anti + anti * inv * t) for t in lifts]
        )
        jac = minus[h][0]
        assert rank(jac.row_join(mixed)) == rank(jac)
    opposite = opposite_cut(w)
    assert rank(opposite) == 10
    left = sp.Matrix.hstack(*kernel(opposite.T))
    right = sp.Matrix.hstack(*kernel(opposite))
    assert rank(left.T * opposite_cut(anti) * right) == 6
    print(
        "Surjective symmetric derivative, mixed obstruction zero, opposite derivative rank 6.",
        flush=True,
    )


def rank_sixteen():
    w, _, second = rational_example()
    assert exchange(w) == second
    first_normal, first_inverse = normal_data(w)
    second_normal, second_inverse = normal_data(second)

    def project_pair(t, u):
        return first_normal(t).col_join(second_normal(u))

    def project(t):
        return project_pair(t, exchange(t))

    data = [tangent(w, h, 0, project) for h in K]
    assert [len(lifts) for _, lifts, _ in data] == [5, 2, 4, 4]
    assert sum(dim for _, _, dim in data) == 61
    for h in K:
        tests = []
        for i, j in combinations_with_replacement(K, 2):
            if i ^ j != h:
                continue
            for k, s in enumerate(data[i][1]):
                for m, t in enumerate(data[j][1]):
                    if i == j and k > m:
                        continue
                    v, u = exchange(s), exchange(t)
                    tests.append(
                        project_pair(
                            s * first_inverse * t + t * first_inverse * s,
                            v * second_inverse * u + u * second_inverse * v,
                        )
                    )
        jac = data[h][0]
        assert rank(jac.row_join(sp.Matrix.hstack(*tests))) == rank(jac)
    print(
        "Both adjacent cuts: all 120 polarized primary obstructions vanish.", flush=True
    )


if __name__ == "__main__":
    rank_ten()
    rank_sixteen()
