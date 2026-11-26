"""Point configurations and linearized projectivities over the split prime 1009."""
from itertools import combinations, combinations_with_replacement, permutations, product

import numpy as np

P = 1009
I = np.eye(3, dtype=np.int64)
FRAME = np.eye(4, dtype=np.int64)
PAIRS = [(1, 2), (2, 0), (0, 1)]


def reduced(matrix):
    matrix = np.array(matrix, dtype=np.int64) % P
    pivots = []
    for col in range(matrix.shape[1]):
        row = len(pivots)
        candidates = np.flatnonzero(matrix[row:, col])
        if not len(candidates):
            continue
        pivot = row + candidates[0]
        matrix[[row, pivot]] = matrix[[pivot, row]]
        matrix[row] = matrix[row] * pow(int(matrix[row, col]), -1, P) % P
        factors = matrix[:, col].copy()
        factors[row] = 0
        matrix = (matrix - factors[:, None] * matrix[row]) % P
        pivots.append(col)
    return matrix, pivots


def inverse(matrix):
    n = len(matrix)
    array, pivots = reduced(np.column_stack((matrix, np.eye(n, dtype=np.int64))))
    assert pivots == list(range(n))
    return array[:, n:]


def normalize(vector):
    vector = np.array(vector, dtype=np.int64) % P
    return vector * pow(int(vector[np.flatnonzero(vector)[0]]), -1, P) % P


def kernel(matrix):
    array, pivots = reduced(matrix)
    assert len(pivots) == 3 and array.shape[1] == 4
    free = next(j for j in range(4) if j not in pivots)
    v = np.zeros(4, dtype=np.int64)
    v[free] = 1
    v[pivots] = -array[:3, free]
    return normalize(v)


def cross(v):
    a, b, c = v
    return np.array([[0, -c, b], [c, 0, -a], [-b, a, 0]])


G = np.array([abs(cross(row)) for row in I])
Q = np.array([-cross(row) for row in I])


def value(coeff, net):
    return np.einsum('i,ijk->jk', coeff, net) % P


def edge(b, a, left):
    sign = 1 if left else -1
    return np.r_[(a[0]*I + sign*value(a[1:], G)) @ b[1:],
                 (b[0]*I - sign*value(b[1:], Q)) @ a[1:]] % P


def configurations():
    root = next(i for i in range(P) if (2*i*i + 1) % P == 0)
    odd = [(1, 0, 0, 0)]
    for zero in range(3):
        for x, y in product((root, -root), repeat=2):
            v = [x, y]
            v.insert(zero, 0)
            odd.append((1, *v))
    odd.extend((0, *v) for v in I)
    odd.extend((0, 1, x, y) for x, y in product((1, -1), repeat=2))
    odd = [normalize(v) for v in odd]
    ends = [[kernel(np.stack([edge(b, e, left) for e in FRAME], axis=1))
             for b in odd] for left in (True, False)]
    for points in ends:
        planes = set()
        for triple in combinations(points, 3):
            normal = kernel(triple)  # Also checks there is no collinear triple.
            planes.add(tuple(i for i, x in enumerate(points) if normal @ x % P == 0))
        assert max(map(len, planes)) == 6
        for i, x in enumerate(points):
            kind = tuple(sum(i in h and len(h) == n for h in planes) for n in (4, 6))
            assert kind == ((12, 6) if np.count_nonzero(x) == 1 else (3, 9))
        assert sum(np.count_nonzero(x) == 1 for x in points) == 4
        assert len({tuple(x) for x in points}) == 20
    return odd, *ends


def bracket(a, b):
    return a @ b - b @ a


def solve_net(dg, df):
    assert not np.any(sum(bracket(u, -b) + bracket(a, v)
                          for a, b, u, v in zip(G, G, dg, df)) % P)
    basis = []
    for k in range(8):
        h = np.zeros((3, 3), dtype=np.int64)
        h.flat[k] = 1
        h[2, 2] = -np.trace(h)
        basis.append(h)
    columns = []
    for pos in range(3):
        for h in basis:
            tiles = [int(pos == i)*bracket(h, Q[j]) + int(pos == j)*bracket(Q[i], h)
                     for i, j in PAIRS]
            columns.append(np.concatenate([x.flatten()[:8] for x in tiles]))
    targets = []
    for i, j in PAIRS:
        target = [[-(bracket(dg[a], -G[c]) + bracket(G[a], df[c]))[i, j]
                   for c in range(3)] for a in range(3)]
        targets.extend(np.array(target).flatten()[:8])
    coords = inverse(np.array(columns).T) @ targets % P
    return [value(row, basis) for row in coords.reshape(3, 8)]


def point_jets(odd, points, left, dg, df, dq):
    even_jets, odd_jets = [], []
    for b, a in zip(odd, points):
        ib, ia = [[i for i in range(4) if i != np.flatnonzero(v)[0]] for v in (b, a)]
        jac = np.stack([edge(e, a, left) for e in FRAME[ib]] +
                       [edge(b, e, left) for e in FRAME[ia]], axis=1)
        q = value(b[1:], dq)
        source = np.r_[value(a[1:], dg if left else df) @ b[1:],
                       (q.T if left else q) @ a[1:]]
        deriv = -inverse(jac) @ source % P
        db, da = np.zeros((2, 4), dtype=np.int64)
        db[ib], da[ia] = deriv[:3], deriv[3:]
        odd_jets.append(db)
        even_jets.append(da)
    return odd_jets, even_jets


def candidates(points, target):
    lookup = {tuple(row): i for i, row in enumerate(target)}
    pivot = next(x for x in points if all(x))
    for axes in permutations(range(4)):
        perm = FRAME[list(axes)]
        reciprocals = np.array([pow(int(x), -1, P) for x in perm @ pivot])
        for y in target:
            if all(y):
                m = (y * reciprocals % P)[:, None] * perm
                m = normalize(m.flatten()).reshape(4, 4)
                labels = [lookup.get(tuple(normalize(m @ x))) for x in points]
                if None not in labels:
                    yield m, labels


def projectivity_equations(points, target, m, labels, dx, dy):
    lhs, rhs = [], []
    for x, vx, j in zip(points, dx, labels):
        y = target[j]
        k = np.flatnonzero(y)[0]
        delta = ((m @ x)[k] * dy[j] - m @ vx) % P
        for h in range(4):
            if h != k:
                e = FRAME[h] - y[h]*FRAME[k]
                lhs.append(np.outer(e, x).flatten())
                rhs.append(e @ delta)
    return np.array(lhs), np.array(rhs) % P


def main():
    odd, a, c = configurations()

    def jets(dg, df):
        dq = solve_net(dg, df)
        left, da = point_jets(odd, a, True, dg, df, dq)
        right, dc = point_jets(odd, c, False, dg, df, dq)
        assert np.array_equal(left, right)
        return da, dc

    dg = np.array([value(row, G) for row in [[3, 2, -1], [1, -3, 4], [-1, 2, 4]]])
    da, dc = jets(dg, -dg)
    for target, dy, survivor in [(a, da, FRAME), (c, dc, np.diag([1, -1, -1, -1]) % P)]:
        maps = list(candidates(a, target))
        assert len(maps) == 96
        kept = []
        for m, labels in maps:
            lhs, rhs = projectivity_equations(a, target, m, labels, da, dy)
            assert len(reduced(lhs)[1]) == 15
            if len(reduced(np.column_stack((lhs, rhs)))[1]) == 15:
                kept.append((m, labels))
        assert len(kept) == 1
        m, labels = kept[0]
        assert np.array_equal(m, survivor)

    # Keep G fixed. First vary F within its span, then in its diagonal quotient.
    def rhs(df):
        dx, dy = jets(np.zeros_like(G), df)
        return projectivity_equations(a, c, m, labels, dx, dy)[1]

    lhs = projectivity_equations(a, c, m, labels, da, dc)[0]
    columns = list(lhs.T)
    for i, j in combinations_with_replacement(range(3), 2):
        k = np.zeros((3, 3), dtype=np.int64)
        k[i, j] = k[j, i] = 1
        columns.append(rhs([value(row, G) for row in k]))
    assert len(reduced(np.array(columns).T)[1]) == 15
    for i in range(3):
        for diag in ([1, 0, -1], [0, 1, -1]):
            df = np.zeros_like(G)
            df[i] = np.diag(diag)
            skew = sum(bracket(x, y) for x, y in zip(G, df))
            for s, t in PAIRS:
                entry = bracket(G[s], G[t])[s, t]
                df[s] -= int(skew[s, t])*int(entry)*G[t]
                assert abs(entry) == 1
            columns.append(rhs(df))
    assert len(reduced(np.array(columns).T)[1]) == 21
    print("PASS: 96 candidates per pair, identity and J survive the mixed tangent; transverse rank 6")


if __name__ == '__main__':
    main()
