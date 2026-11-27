"""Opposite determinant and four monic top tensors on the extra simple component."""
from itertools import product

import sympy as sp

from verify_split_inverse_windows import bell
from verify_split_simple_surface import EXTRA, D0
from verify_split_spectral_line import multiply, reconstruct


def check():
    t = sp.Symbol('t')
    tensor = reconstruct(bell((sp.Integer(1), -sp.I*t, sp.I, t)), multiply(EXTRA, D0))
    words = list(product(range(4), repeat=4))
    orders = [(0, 1, 2, 3), (1, 0, 2, 3), (0, 1, 3, 2), (1, 0, 3, 2)]
    volumes, pairings = [], []
    for i in range(4):
        def value(q):
            return tensor.get(q[-i:] + q[:-i] if i else q, 0)

        def desc(q, start):
            labels = [orders[(start+j) % 4][a] for j, a in enumerate(q)]
            return all(a > b for a, b in zip(labels, labels[1:]))

        volume = value(next(q for q in words if desc(q, i)))
        volumes.append(sp.factor(volume))
        triples = [q for q in product(range(4), repeat=3) if desc(q, i+1)]
        pairing = sp.Matrix([[value((a,) + q)/volume for q in triples] for a in range(4)])
        pairings.append(sp.factor(pairing.det()))
        pairs = list(product(range(4), repeat=2))
        ordered = sorted(pairs, key=lambda q: tuple(orders[(i+j)%4][a] for j, a in enumerate(q)),
                         reverse=True)
        relations = sp.Matrix([[value(q+p) for q in ordered] for p in pairs])
        # Prioritize the six reducible pairs to check pure binary and all mixed rewrites.
        descending = [q for q in ordered if desc(q, i)]
        ascending = [q for q in ordered if not desc(q, i)]
        columns = [ordered.index(q) for q in descending+ascending]
        reduced, pivots = relations[:, columns].rref()
        assert pivots == tuple(range(6))
        for j, q in enumerate(descending):
            for k, p in enumerate(ascending):
                if reduced[j, 6+k]:
                    assert tuple(orders[(i+b)%4][a] for b,a in enumerate(p)) < tuple(
                        orders[(i+b)%4][a] for b,a in enumerate(q))
    assert sp.factor(sp.prod(pairings)) == 1
    assert all(sp.cancel(a-b) == 0 for a,b in zip(volumes, [-(t-1)*(t+1)/(2*t), 1, 1, 1]))
    opposite = sp.Matrix(16, 16, [tensor.get((a, c, b, d), 0) for a, b, c, d in words])
    determinant = sp.factor(opposite.det())
    assert sp.cancel(determinant -
                     (t*t+1)**4*(t*t-2*t-1)**2*(t*t+2*t-1)**2/(256*t**8)) == 0


if __name__ == '__main__':
    check()
    print('Four PBW top volumes, monic pairing product and opposite determinant verified.')
