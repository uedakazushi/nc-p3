"""Linear parameters and shift on marked indexed Ore extensions of a quadratic AS3 base."""

import sympy as sp


def span(matrix):
    basis = matrix.columnspace()
    return sp.Matrix.hstack(*basis) if basis else sp.zeros(matrix.rows, 0)


def kernel(matrix):
    basis = matrix.nullspace()
    return sp.Matrix.hstack(*basis) if basis else sp.zeros(matrix.cols, 0)


def quotient(relations):
    """Quotient projection and representatives in the standard tensor basis."""
    n = relations.rows
    basis = span(relations.row_join(sp.eye(n)))
    k = relations.rank()
    return basis.inv()[k:, :], basis[:, k:]


def intersection(left, right):
    return span(left * kernel(left.row_join(-right))[: left.cols, :])


def parameters(relations):
    identity = sp.eye(3)
    cubic = sp.kronecker_product(relations, identity).row_join(
        sp.kronecker_product(identity, relations)
    )
    q2, lift = quotient(relations)
    q3, _ = quotient(cubic)
    assert q2.shape == (6, 9) and q3.shape == (10, 27)
    left, right = [], []
    for j in range(3):
        for a in range(6):
            op = lift[:, a] * identity[j, :]
            left.append(
                (q3 * sp.kronecker_product(op, identity) * relations).reshape(30, 1)
            )
            right.append(
                (q3 * sp.kronecker_product(identity, op) * relations).reshape(30, 1)
            )
    left, right = sp.Matrix.hstack(*left), sp.Matrix.hstack(*right)
    assert left.rank() == right.rank() == 15
    translations = []
    for j in range(3):
        for op in [
            sp.kronecker_product(identity[:, j], identity),
            sp.kronecker_product(identity, identity[:, j]),
        ]:
            translations.append((q2 * op).T.reshape(18, 1))
    shifts = span(sp.Matrix.hstack(*translations))
    e = sp.eye(18)
    dimensions = [18]
    while True:
        pairs = kernel((left * e).row_join(right * e))
        first, second = e * pairs[: e.cols, :], e * pairs[e.cols :, :]
        new = intersection(first, second)
        if new.cols == e.cols:
            break
        e = new
        dimensions.append(e.cols)
    basis = span(shifts.row_join(e))
    assert basis.cols == e.cols
    pivots = basis.T.rref()[1]
    proj = basis[list(pivots), :].inv()[shifts.cols :, :]
    a, b = [proj * m[list(pivots), :] for m in (first, second)]
    indices = a.rref()[1]
    shift = b[:, list(indices)] * a[:, list(indices)].inv()
    assert shift * a == b and shift.rank() == shift.rows
    return dimensions, shifts.cols, shift
