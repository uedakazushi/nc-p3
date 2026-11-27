"""Full thin rank-one cross blocks in the bisymmetric Klein character type."""
from itertools import combinations, product

import sympy as s
from flint import fmpz_mat

from verify_negative_component_quotient import K, relabel, thin, words
from verify_parallel_invertible_rectangles import smith

X, Y, W, V = parameters = s.symbols("X Y W V")


def coefficients():
    full = s.Matrix([[X, X, 1, 1], [Y, Y, 1, 1], [1, 1, W, W], [1, 1, V, V]])
    blocks = full, full.T, s.ones(4), s.ones(4)
    return {
        q: blocks[q[0] ^ q[1]][q[0], q[2]]
        * (-1)**sum((q[i] >> 1) * (q[j] & 1) for i, j in combinations(K, 2))
        for q in words
    }


def index_row(q):
    return [int(q[p] == i) for p in K for i in K]


def check_gauge():
    rows = [index_row(q) for q in thin]
    base = fmpz_mat(rows).hnf().tolist()[:12]
    top, bottom = [index_row((i,) * 4) for i in (0, 2)]
    differences = [[a + b for a, b in zip(top, bottom)]]
    for q in set(words) - set(thin):
        ref = top if max(q) < 2 else bottom
        differences.append([a - b for a, b in zip(index_row(q), ref)])
    for diff in differences:
        assert fmpz_mat(rows + [diff]).hnf().tolist()[:12] == base
    diag = fmpz_mat(rows + [top]).snf()
    assert [diag[i, i] for i in range(16)] == [1]*11 + [2, 2, 0, 0, 0]


def check_permutations(F):
    E = s.Matrix([index_row(q) for q in thin])
    H, L, R = smith(E)
    rank = E.rank()
    picks = [next(q for q in words if F[q] in (x, -x)) for x in parameters]
    images = []
    for shear, i in product(range(2), K):
        for j, k in product((i & 2, (i & 2) ^ 1), repeat=2):
            G = {q: F[relabel(q, 0, shear, (i, j, k, i ^ j ^ k))] for q in words}
            rhs = L * s.Matrix([2 * int(G[q] != F[q]) for q in thin])
            assert all(v % 4 == 0 for v in rhs[rank:])
            factors = R * s.Matrix([rhs[t]/H[t, t] if t < rank else 0 for t in range(16)])
            assert all(v.is_integer for v in factors)
            for q in words:
                G[q] *= s.I**(int(sum(factors[4*p + n] for p, n in enumerate(q))) % 4)
            assert all(G[q] == F[q] for q in thin)
            sub = {x: s.cancel(x*G[q]/F[q]) for x, q in zip(parameters, picks)}
            defects = {s.expand(G[q] - F[q].xreplace(sub)) for q in words} - {0}
            if defects:
                assert defects <= {X-Y, Y-X, W-V, V-W}
                assert defects & {X-Y, Y-X} and defects & {W-V, V-W}
            else:
                sign = -1 if sub[X].could_extract_minus_sign() else 1
                images.append(tuple(sign*v for v in sub.values()))
    target = {left + right for lo, hi in [((X, Y), (W, V)), ((W, V), (X, Y))]
              for left in [lo, lo[::-1]] for right in [hi, hi[::-1]]}
    assert len(images) == 8 and set(images) == target


def check_ranks(F):
    opposite = [s.Matrix(4, 4, lambda i, j: F[i, j, i ^ h, j ^ h]) for h in K]
    assert all(s.expand(m.det() - d) == 0 for m, d in
               zip(opposite, [-4*(V+W)*(X+Y)]*2 + [-16]*2))
    for m in opposite[:2]:
        for sub, rank in [({Y: -X}, 3), ({V: -W}, 3), ({Y: -X, V: -W}, 2)]:
            mat = m.subs(sub)
            minors = [s.expand(mat.extract(I, J).det()) for I in combinations(K, rank)
                      for J in combinations(K, rank)]
            # The lower bound is uniform after excluding the other vanishing sum.
            expected = (V+W if Y in sub else X+Y) if rank == 3 else 1
            assert s.factor(s.gcd_list(minors)/expected) in (1, -1, 2, -2, 4, -4)
            assert mat.rank() == rank
    A, B, C, D, r, t = s.symbols("A B C D r t")
    web = s.Matrix([[A, 0, r, t], [0, B, -t, -r], [r, -t, C, 0], [t, -r, 0, D]])
    minors = {s.factor(web.extract(I, J).det()) for I in combinations(K, 3)
              for J in combinations(K, 3)}
    for e in (r*t*(A+B), r*t*(C+D)):
        assert e in minors or -e in minors


if __name__ == "__main__":
    F = coefficients()
    check_gauge()
    check_permutations(F)
    check_ranks(F)
    print("PASS: all singular-corner sign corrections, thick weights, ranks and quadric minors")
