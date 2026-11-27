"""Exact cut, boundary, congruence and stabilizer calculations for the single-standard locus."""
from itertools import product

import sympy as sp

from verify_sl2_window_profiles import EPS, K, repeated_tensor

WORDS = list(product(range(2), repeat=4))


def ordinary_equations(boundaries):
    rows = []
    for i in (0, 1):
        annihilators = [boundaries[j].reshape(1, 4).nullspace() for j in (i, i + 2)]
        for left, right in product(*annihilators):
            rows.append([left[2 * w[i] + w[i + 1]] * right[2 * w[(i + 2) % 4] + w[(i + 3) % 4]]
                         for w in WORDS])
    return sp.Matrix(rows)


def check_singular():
    # A full projection is recorded separately at each end of each mixed plane.
    possible = set()
    for flags in product((False, True), repeat=8):
        left, right = flags[:4], flags[4:]
        if not all(left[i] or right[i] for i in K):
            continue
        if not all(left[i] or right[(i - 1) % 4] for i in K):
            continue
        opposite = []
        for i in (0, 1):
            a = left[i] and left[i + 2]
            b = right[(i - 1) % 4] and right[(i + 1) % 4]
            if a != b:
                break
            opposite.append(a)
        if len(opposite) != 2:
            continue
        regular_edges = tuple(right[(i - 1) % 4] and left[(i + 1) % 4] for i in K)
        possible.add((sum(opposite), min(regular_edges[i:] + regular_edges[:i] for i in K)))
    assert possible == {(1, (False, False, True, True)), (2, (True,) * 4)}
    for first in (sp.diag(1, 0), sp.Matrix([[0, 1], [0, 0]])):
        maps = (first, sp.diag(1, 0), sp.eye(2), sp.eye(2))
        assert ordinary_equations([EPS * m.T for m in maps]).rank() == 16


def check_ordinary():
    t = sp.symbols('t', nonzero=True)
    def equations(monodromy):
        return ordinary_equations([EPS] * 3 + [EPS * monodromy.T])
    mat = equations(sp.diag(1, t))
    rows = (0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16)
    assert sp.factor(mat[list(rows), :].det() + (t - 1)**3 / t**3) == 0
    assert equations(sp.Matrix([[1, 1], [0, 1]])).rank() == 16
    basis = sp.Matrix.hstack(*equations(sp.eye(2)).nullspace())
    assert basis.shape == (16, 3)
    assert sp.Matrix(16, 3, lambda i, j: basis[i % 8 * 2 + i // 8, j]) == -basis
    a, b, c, d = sp.symbols('a b c d')
    g = sp.Matrix([[a, b], [c, d]])
    sym = sp.Matrix([[a*a, a*b, b*b], [2*a*c, a*d+b*c, 2*b*d], [c*c, c*d, d*d]])
    assert (sp.kronecker_product(g, g, g, g) * basis - g.det() * basis * sym).applyfunc(sp.expand).is_zero_matrix
    p = g
    monodromy = -p * p.inv().T
    q = -EPS * p.inv().T * EPS
    assert (q + EPS * p.inv() * (EPS * monodromy.T).T).applyfunc(sp.cancel).is_zero_matrix
    assert (q - p / p.det()).applyfunc(sp.cancel).is_zero_matrix


def projected_dimension(w):
    # The zero-action equation already permits all projective stabilizer directions
    # because scalar summands can cancel an arbitrary multiple of w.
    columns = [sp.Matrix([w[word[:i] + (b,) + word[i+1:]] if word[i] == a else 0
                          for word in product(K, repeat=4)]) for i, a, b in product(K, K, K)]
    lifts = sp.Matrix.hstack(*columns).nullspace()
    first = [c[:16, :].reshape(4, 4) for c in lifts]
    return sp.Matrix.hstack(*[(m - sp.eye(4)*m.trace()/4).reshape(16, 1) for m in first]).rank()


def check_other_profiles():
    def adjoint(q):
        zeros = [i for i, v in enumerate(q) if not v]
        if len(zeros) == 1:
            i = zeros[0]
            return sp.LeviCivita(*(q[(i+j) % 4] for j in (1, 2, 3)))
        if zeros == [0, 2]:
            return int(q[1] == q[3])
        if zeros == [1, 3]:
            return -int(q[0] == q[2])
        return 0
    assert projected_dimension({q: adjoint(q) for q in product(K, repeat=4)}) == 3
    maps = (sp.eye(2),)*3 + (sp.Matrix([[1, 1], [0, 1]]),)
    assert projected_dimension(repeated_tensor(maps, 1, -1)) == 4


if __name__ == '__main__':
    check_singular()
    check_ordinary()
    check_other_profiles()
    print('Singular boundaries excluded. Scalar monodromy has a binary quadratic slice; all other slices vanish.')
