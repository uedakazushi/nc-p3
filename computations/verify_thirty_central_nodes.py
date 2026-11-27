"""Certify nodality on the fifteen-dimensional component using an integral local chart."""
import itertools
import math
import random
from functools import cache
from multiprocessing import get_context
import numpy as np
import sympy as sp
from flint import nmod_mat, nmod_mpoly_ctx
import verify_quartic_center_order as original
from verify_klein_elliptic_completion import rational_example, four_relations
P = 1009
DIGITS = 5
MOD = P ** DIGITS
assert MOD < 2 ** 63


def array(m):
    return np.array(m.tolist(), dtype=np.int64)


def mat(a, m):
    return nmod_mat(a.tolist() if isinstance(a, np.ndarray) else a, m)


def piv(m):
    red, rank = m.rref()
    arr = array(red)
    return [int(np.flatnonzero(row)[0]) for row in arr[:rank]]


def inv(m, mod):
    """Newton lift the inverse of a matrix whose determinant is a unit."""
    a = array(m)
    x = mat(array(mat(a % P, P).inv()), mod)
    ident = mat(np.eye(len(a), dtype=np.int64) * 2, mod)
    for _ in range((DIGITS - 1).bit_length()):
        x = x * (ident - m * x)
    return x


def unit_reduction(m):
    if not m.nrows():
        return (m, 0)
    tab = array(m)
    base = tab % P
    cols = piv(mat(base, P))
    if not cols:
        return m, 0
    rows = piv(mat(base[:, cols].T, P))
    r = inv(mat(tab[np.ix_(rows, cols)], MOD), MOD) * mat(tab[rows], MOD)
    return (r, len(cols))


def normal_projection(m):
    rows, r = unit_reduction(m)
    reduced = array(rows)
    pivots = piv(mat(array(m) % P, P))
    free = [j for j in range(m.ncols()) if j not in pivots]
    result = np.zeros((len(free), m.ncols()), dtype=np.int64)
    result[:, free] = np.eye(len(free), dtype=np.int64)
    result[:, pivots] = -reduced[:r, free].T % MOD
    q = mat(result, MOD)
    assert not any((m * q.transpose()).entries())
    return (q, free)


def residue(x):
    y = sp.Rational(x)
    return int(y.p) * pow(int(y.q), -1, MOD) % MOD


def lift_window():
    """Move the rational window in an integral chart of the two rank-six loci."""
    w, _, rev = rational_example()
    base = np.array([[residue(v) for v in row] for row in w.tolist()], dtype=object).reshape((4,) * 4)
    groups = []
    jac = []
    rng = random.Random(20883)
    for axes in ((0, 1, 2, 3), (1, 2, 3, 0)):
        a = np.array(base.transpose(axes).reshape(16, 16) % P, dtype=np.int64)
        ii = piv(mat(a.T, P))
        jj = piv(mat(a, P))
        uu = [i for i in range(16) if i not in ii]
        vv = [j for j in range(16) if j not in jj]
        mi = array(mat(a[np.ix_(ii, jj)], P).inv())
        left = np.zeros((10, 16), dtype=np.int64)
        right = left.copy()
        left[:, uu] = np.eye(10, dtype=np.int64)
        right[:, vv] = np.eye(10, dtype=np.int64)
        left[:, ii] = -(a[np.ix_(uu, jj)] @ mi) % P
        right[:, jj] = -(mi @ a[np.ix_(ii, vv)]).T % P
        jac.extend((np.outer(l, r).reshape((4,) * 4).transpose(np.argsort(axes)).ravel() % P for l in left for r in right))
        groups.append((axes, ii, jj, uu, vv))
    jac = np.array(jac)
    rows = piv(mat(jac.T, P))
    cols = piv(mat(jac, P))
    assert len(rows) == len(cols) == 180
    solve = array(mat(jac[np.ix_(rows, cols)], P).inv())
    direction = np.array([rng.randrange(P) for _ in range(256)], dtype=np.int64)
    direction[cols] -= solve @ (jac[rows] @ direction % P) % P
    current = base + P * direction.reshape((4,) * 4)
    for n in range(2, DIGITS + 1):
        residual = []
        for axes, ii, jj, uu, vv in groups:
            a = np.array(current.transpose(axes).reshape(16, 16) % MOD, dtype=np.int64)
            c = mat(a[np.ix_(uu, vv)], MOD) - mat(a[np.ix_(uu, jj)], MOD) * inv(mat(a[np.ix_(ii, jj)], MOD), MOD) * mat(a[np.ix_(ii, vv)], MOD)
            residual.extend(map(int, c.entries()))
        r = np.array(residual, dtype=np.int64)
        assert not (r % P ** n).any()
        if n < DIGITS:
            delta = np.zeros(256, dtype=np.int64)
            delta[cols] = -solve @ (r[rows] // P ** n % P) % P
            current += P ** n * delta.astype(object).reshape((4,) * 4)
    rules = []
    for i in range(4):
        a = np.array(current.transpose([(i + j) % 4 for j in range(4)]).reshape(16, 16) % MOD, dtype=np.int64)
        basis = a[:, piv(mat(a % P, P))].T
        rules.append([({divmod(j, 4): int(v) for j, v in enumerate(row) if v}, None) for row in basis])
    return (rules, original.integral_quartics(four_relations(w, rev), P)[1])


def lift_centers(rules, old):
    """Compute the path bases and lift the five commuting quartics in them."""
    bases = {(i, 0): [()] for i in range(4)}
    maps = {}

    @cache
    def image(i, word):
        if not word:
            return mat([[1]], MOD)
        return maps[i, len(word), word[0]] * image((i + 1) % 4, word[1:])
    for n in range(1, 6):
        for i in range(4):
            suffix = bases[(i + 1) % 4, n - 1]
            words = [(j,) + v for j in range(4) for v in suffix]
            rows = []
            if n >= 2:
                for rule, _ in rules[i]:
                    for tail in bases[(i + 2) % 4, n - 2]:
                        a = [nmod_mat(len(suffix), 1, MOD) for _ in range(4)]
                        for (j, k), c in rule.items():
                            a[j] += c * image((i + 1) % 4, (k,) + tail)
                        rows.append([int(v) for block in a for v in block.entries()])
            if rows:
                q, free = normal_projection(mat(rows, MOD))
            else:
                q = mat(np.eye(4, dtype=np.int64), MOD)
                free = list(range(4))
            bases[i, n] = [words[g] for g in free]
            for j in range(4):
                maps[i, n, j] = mat(array(q)[:, len(suffix) * j:len(suffix) * (j + 1)], MOD)
    eq = []
    for i, j in itertools.product(range(4), repeat=2):
        k = (i + 1) % 4
        left = np.array([[int(v) for v in image(i, w + (j,)).entries()] for w in bases[i, 4]], dtype=np.int64).T
        right = np.array([[int(v) for v in image(i, (j,) + w).entries()] for w in bases[k, 4]], dtype=np.int64).T
        rows = np.zeros((56, 140), dtype=np.int64)
        rows[:, 35 * i:35 * (i + 1)] = left
        rows[:, 35 * k:35 * (k + 1)] = -right
        eq.extend(rows)
    ker, free = normal_projection(mat(np.array(eq) % MOD, MOD))
    assert len(free) == 5
    out = []
    for z in old:
        v = np.concatenate([sum((int(c) * array(image(i, w))[:, 0].astype(object) for w, c in z[i].items())) % P for i in range(4)])
        coeff = [int(c) for c in (mat([list(map(int, v[free]))], MOD) * ker).entries()]
        assert list(v) == [i % P for i in coeff]
        out.append([{w: c for w, c in zip(bases[i, 4], coeff[35 * i:35 * (i + 1)]) if c} for i in range(4)])
    return out


class LiftedFiber(original.CentralFiber):

    def pivot_columns(self, constraint):
        reduced, rank = unit_reduction(constraint)
        columns = [next((j for j in range(reduced.ncols()) if int(reduced[i, j]) % P)) for i in range(rank)]
        return (reduced, columns)

    def invert(self, matrix):
        return inv(matrix, MOD)


def sample(args):
    rules, quartics, point = args
    f = LiftedFiber(rules, quartics, (1, *point), MOD)
    a = f.fifth_action(quartics[4])
    # The regular characteristic polynomial is the sixteenth power of the central relation.
    traces = []
    acc = a
    for i in range(4):
        traces.append(sum((int(acc[j, j]) for j in range(64))) * pow(16, -1, MOD) % MOD)
        acc = acc * a
    coef = [1]
    for i in range(1, 5):
        coef.append(-sum((coef[j] * traces[i - j - 1] for j in range(i))) * pow(i, -1, MOD) % MOD)
    f.prefix.cache_clear()
    return coef


def central_form():
    rules, old = lift_window()
    quartics = lift_centers(rules, old)
    pts = [pt for pt in itertools.product(range(5), repeat=3) if sum(pt) <= 4]
    with get_context('fork').Pool(4) as pool:
        vals = pool.map(sample, [(rules, quartics, pt) for pt in pts])
    evaluation = [[math.prod(x**i for x, i in zip(pt, e)) for e in pts] for pt in pts]
    coeff = inv(mat(evaluation, MOD), MOD) * mat(vals, MOD)
    a, b, c, d, z = sp.symbols('a b c d z')
    result = 0
    for j, e in enumerate(pts):
        for i in range(5):
            q = int(coeff[j, i])
            assert sum(e) <= i or not q
            if q:
                result += q * z ** (4 - i) * a ** (i - sum(e)) * sp.prod((v ** k for v, k in zip((b, c, d), e)))
    return result


def lower_rank(a, precision):
    """Count independent pivots by valuation, retaining only known digits."""
    tally = []
    while precision and min(a.shape):
        modulus = P ** precision
        cols = piv(mat(a % P, P))
        tally.append(len(cols))
        if not cols:
            a = a // P
            precision -= 1
            continue
        othercols = [i for i in range(a.shape[1]) if i not in cols]
        if not othercols:
            break
        rows = piv(mat((a[:, cols] % P).T, P))
        otherrows = [i for i in range(a.shape[0]) if i not in rows]
        if not otherrows:
            break
        a11 = mat(a[np.ix_(rows, cols)], modulus)
        b = mat(a[np.ix_(rows, othercols)], modulus)
        reduced_inv = mat(a[np.ix_(rows, cols)] % P, P).inv()
        x = array(reduced_inv * mat(array(b) % P, P))
        for n in range(1, precision):
            residual = array(b - a11 * mat(x, modulus))
            assert not (residual % P ** n).any()
            x += P ** n * array(reduced_inv * mat(residual // P ** n % P, P))
        second = mat(a[np.ix_(otherrows, othercols)], modulus) - mat(a[np.ix_(otherrows, cols)], modulus) * mat(x, modulus)
        a = array(second)
        assert not (a % P).any()
        a //= P
        precision -= 1
    return tally


def certify(form):
    vars = sp.symbols('a b c d z')
    ctx = nmod_mpoly_ctx.get(tuple(map(str, vars)), MOD)
    f = ctx.from_dict({k: int(v) % MOD for k, v in sp.Poly(form, *vars).terms()})
    jac = [f.derivative(i) for i in range(5)]
    hess = [[q.derivative(i) for i in range(5)] for q in jac]
    cof = []
    for i, j in itertools.combinations_with_replacement(range(5), 2):
        rows = [u for u in range(5) if u != i]
        cols = [u for u in range(5) if u != j]
        det = ctx.constant(0)
        for perm in itertools.permutations(range(4)):
            sign = (-1) ** sum((perm[i] > perm[j] for i, j in itertools.combinations(range(4), 2)))
            det += sign * math.prod((hess[r][cols[k]] for r, k in zip(rows, perm)))
        cof.append(det)

    def mons(n):
        return list((t for t in itertools.product(range(n + 1), repeat=5) if sum(t) == n))

    def matrix_degree(n, with_cof):
        target = {m: i for i, m in enumerate(mons(n))}
        rows = []
        for pol, weight in [(p, 3) for p in jac] + ([(q, 8) for q in cof] if with_cof else []):
            terms = list(pol.to_dict().items())
            for h in mons(n - weight):
                row = np.zeros(len(target), dtype=np.int64)
                for k, v in terms:
                    row[target[tuple((i + j for i, j in zip(k, h)))]] = int(v) % MOD
                rows.append(row)
        return np.array(rows)
    for with_cofactors, expected in [(False, [637, 0, 48, 0, 0]), (True, [643, 0, 56, 0, 16])]:
        matrix = matrix_degree(9, with_cofactors)
        ranks = lower_rank(matrix, DIGITS)
        assert ranks == expected
        print('Degree 9, Hessian cofactors included:', with_cofactors, 'pivot counts:', ranks, flush=True)


if __name__ == '__main__':
    certify(central_form())
