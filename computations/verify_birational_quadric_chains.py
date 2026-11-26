"""Check the degree-three overlap for switched rulings of both quadric types."""
from itertools import product
from random import Random

from sympy import Matrix


def sections(degrees):
    return [(i, k) for i, d in enumerate(degrees) for k in range(d + 1)]


def check(left, right):
    rng = Random(1729)
    v, w = sections(left), sections(right)
    a = sections([d + 1 for d in left])
    b = sections([d + 1 for d in right])
    basis = list(product(a, b))
    p = {(i, k, t): rng.randint(1, 9) for (i, k), t in product(v, range(2))}
    q = {(s, j, l): rng.randint(1, 9) for s, (j, l) in product(range(2), w)}
    pcols, qcols = [], []
    for s, (j, l) in product(range(2), w):
        pcols.append([p.get((i, k-s, h-l), 0) if j == u else 0
                      for (i, k), (u, h) in basis])
    for (i, k), t in product(v, range(2)):
        qcols.append([q.get((h-k, j, l-t), 0) if i == u else 0
                      for (u, h), (j, l) in basis])
    pm, qm = Matrix(pcols).T, Matrix(qcols).T
    pq = pm * Matrix(list(q.values()))
    assert pq == qm * Matrix(list(p.values())) and any(pq)
    assert len(basis) == 36 and pm.rank() == qm.rank() == 8
    assert pm.row_join(qm).rank() == 15


for left, right in product(((1, 1), (2, 0)), repeat=2):
    check(left, right)
print("PASS: both quadric types in both orders give a common product and cubic dimension 21")
