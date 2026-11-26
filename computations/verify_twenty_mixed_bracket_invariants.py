"""Rational mixed-bracket frames and re-marking on the 20-sheet cover."""
import numpy as np
import sympy as sp
from verify_mixed_bracket_reconstruction import (
    INDICES, P, brackets, kernel_over_dvr, lift, operator, trace_matrix, units,
)
from verify_semisimple_period_obstruction import configurations, inverse, reduced

DIGITS = 14
MOD = P**DIGITS
IDENTITY = np.eye(3, dtype=object)
BASIS = units(MOD)
DUAL = np.array([IDENTITY[:, j:j+1]@IDENTITY[i:i+1, :]-int(i == j)*np.diag([0, 0, 1])
                 for i, j in INDICES], dtype=object)


def tensor(data):
    g, f, q, p = data
    w = np.zeros((4, 4, 4, 4), dtype=object)
    for i in range(3):
        w[0, i+1, 0, i+1] = w[i+1, 0, i+1, 0] = 1
        w[i+1, 1:, 0, 1:], w[0, 1:, i+1, 1:] = g[i], f[i]
        w[1:, i+1, 1:, 0], w[1:, 0, 1:, i+1] = q[i], p[i]
    w[1:, 1:, 1:, 1:] = ((g[:, None]@f[None, :]).transpose(0, 2, 1, 3)
                           +(q[:, None]@p[None, :]).transpose(2, 0, 3, 1))
    return w % MOD


def adjacent_evaluations(w, points):
    return np.concatenate([np.kron(points[i], points[(i+1) % 4])
                           @w.transpose([(i+j) % 4 for j in range(4)]).reshape(16, 16)
                           for i in range(4)])


def apply_frame(w, frame, axis):
    return np.moveaxis(np.tensordot(frame, w, axes=(1, axis)), 0, axis) % MOD


def lift_mark(w, initial):
    points = np.array(initial, dtype=object)
    fixed = [next(j for j in range(4) if points[i, j]) for i in range(4)]
    free = [(i, j) for i, j in np.ndindex(4, 4) if j != fixed[i]]
    start = adjacent_evaluations(w, points)
    cols = []
    for i, j in free:
        shifted = points.copy()
        shifted[i, j] += 1
        cols.append((adjacent_evaluations(w, shifted)-start) % P)
    jac = np.array(cols, dtype=np.int64).T
    rows = reduced(jac.T)[1]
    assert len(rows) == 12
    solve = inverse(jac[rows])
    for n in range(1, DIGITS):
        residual = adjacent_evaluations(w, points)
        assert not np.any(residual % P**n)
        delta = -solve@(residual[rows]//P**n % P).astype(np.int64) % P
        for h, (i, j) in zip(delta, free):
            points[i, j] += int(h)*P**n
    assert not np.any(adjacent_evaluations(w, points) % MOD)
    t = w
    for i in range(4):
        frame = np.vstack([points[i], np.delete(np.eye(4, dtype=object), fixed[i], axis=0)])
        t = apply_frame(t, frame, i)
    for i, gram in [(2, t[1:, 0, 1:, 0]), (3, t[0, 1:, 0, 1:])]:
        m = sp.Matrix(gram)
        frame = np.eye(4, dtype=object)
        frame[1:, 1:] = np.array(m.adjugate().T, dtype=object)*pow(int(m.det()), -1, MOD) % MOD
        t = apply_frame(t, frame, i)
    q = t[1:, 1:, 1:, 0].transpose(1, 0, 2)
    p = t[1:, 0, 1:, 1:].transpose(2, 0, 1)
    c = brackets(q, p) % MOD
    reverse = brackets(t[1:, 1:, 0, 1:], t[0, 1:, 1:, 1:].transpose(1, 0, 2))
    assert not np.any((c.transpose(2, 3, 0, 1)+reverse) % MOD)
    return c


def covariant(operator_matrix):
    image = (operator_matrix.T@BASIS.reshape(8, 9)).reshape(8, 3, 3) % MOD
    a = sum(x@y+y@x for x, y in zip(image, DUAL)) % MOD
    return (a-pow(3, -1, MOD)*a.trace()*IDENTITY) % MOD


def check_reconstruction(c):
    left, n = kernel_over_dvr(trace_matrix(c, BASIS), DIGITS)
    right, m = kernel_over_dvr(trace_matrix(c.swapaxes(0, 1), BASIS), DIGITS)
    def matrices(rows):
        return (rows.reshape(6, 3, 8)@BASIS.reshape(8, 9)).reshape(6, 3, 3, 3) % MOD
    cols = [brackets(a, b).reshape(-1) for a in matrices(left) for b in matrices(right)]
    _, remaining = kernel_over_dvr(cols, min(n, m))
    assert remaining >= DIGITS-4


def scores(c):
    transposed = c.transpose(2, 3, 0, 1)
    check_reconstruction(c)
    check_reconstruction(transposed)
    first, second = operator(c, BASIS), operator(transposed, BASIS)
    result, generators = [], []
    for square in (second@first % MOD, first@second % MOD):
        norm = square*pow(int(square.trace()), -1, MOD) % MOD
        a, b = covariant(norm), covariant(norm@norm % MOD)
        words = [IDENTITY, a, b, a@a, a@b, b@a, b@b, a@a@b, a@b@b]
        _, remaining = kernel_over_dvr(np.array(words).reshape(9, 9), DIGITS)
        assert remaining >= DIGITS-4
        determinant = int(sp.Matrix(a@b-b@a).det()) % MOD
        assert determinant % P**10 != 0
        generators.append(a)
        result.append(int((a@a).trace()) % MOD)
    image = (transposed.reshape(9, 9)@generators[0].T.flatten()).reshape(3, 3)
    multiplier = int((generators[1]@image).trace()) % MOD
    assert multiplier % P**3 == 0 and multiplier % P**4 != 0
    return result


def main():
    assert not np.any((BASIS.reshape(8, 9)@DUAL.transpose(0, 2, 1).reshape(8, 9).T
                       -np.eye(8, dtype=object)) % MOD)
    w = tensor(lift())
    odd, first, third = configurations()
    values = []
    for mark in zip(first, odd, third, odd):
        values.extend(scores(lift_mark(w, mark)))
    assert len({h % P**3 for h in values}) == 40
    print('PASS: all 20 lifted marks, both inverses at each, 40 separated rational values,'
          ' nonzero oriented determinants and 40 nine-word frames')


if __name__ == '__main__':
    main()
