"""The 71-node mixed-bracket bicubic and its equisingular tangent space."""
from itertools import combinations, combinations_with_replacement, product
import os
import subprocess

import numpy as np
import sympy as s
from flint import nmod_mat

from verify_involution_quartet_component import K, kernel, rows
from verify_klein_elliptic_completion import rational_example


PRIME = 2657


def marked_triples():
    p = PRIME
    W, _, _ = rational_example()
    raw = np.array([int(v.p)*pow(int(v.q), -1, p)%p for v in W], dtype=np.int64).reshape([4]*4)
    start = np.array([1, 1, 0, 0], dtype=np.int64)
    point = start
    t = raw
    for i in K:
        pivot = next(j for j in K if point[j])
        frame = np.vstack([point, np.delete(np.eye(4, dtype=np.int64), pivot, axis=0)])
        t = np.moveaxis(np.tensordot(frame, t, axes=(1, i))%p, 0, i)
        local = raw.transpose(tuple((i+j)%4 for j in K))
        nxt = kernel(rows(np.tensordot(point, local, axes=(0, 0)).reshape(4, 16).T, p))
        assert len(nxt) == 1
        point = nxt[0]*pow(int(next(v for v in nxt[0] if v)), -1, p)%p
    assert np.array_equal(point, start)
    for pos, gram in ((2,t[1:, 0, 1:, 0]), (3,t[0, 1:, 0, 1:])):
        frame = np.eye(4, dtype=np.int64)
        frame[1:, 1:] = np.array(rows(gram, p).inv().transpose().tolist(), dtype=np.int64)
        t = np.moveaxis(np.tensordot(frame, t, axes=(1, pos))%p, 0, pos)
    def center(a):
        return (a-np.array([m.trace()*pow(3,-1,p)*np.eye(3,dtype=np.int64) for m in a]))%p
    return np.array(list(map(center, (t[1:,1:,0,1:], t[0,1:,1:,1:].transpose(1,0,2),
                                      t[1:,1:,1:,0].transpose(1,0,2), t[1:,0,1:,1:].transpose(2,0,1)))))


def bracket(A, B):
    return (np.einsum('iab,jbc->ijac', A, B)-np.einsum('iab,jbc->jiac', B, A))%PRIME


def checked_chart(matrix, l, r, i, j, interpolate):
    variables = [z for z in (*l, *r) if z not in (l[i], r[j])]
    fixed = {l[i]: 1, r[j]: 1}
    M = matrix.subs(fixed)
    polynomial = s.Poly(M.det(), *variables, modulus=PRIME).as_expr()
    fmt = lambda f: str(f).replace('**', '^')
    minors = [M.extract(I, J).det() for I, J in product(combinations(range(3),2), repeat=2)]
    # Keep all four derivatives on the affine chart before imposing earlier zero coordinates.
    earlier = [*l[:i], *r[:j]]
    ring = ','.join(map(str, variables))
    program = (f'ring R={PRIME},({ring}),dp;short=0;poly F={fmt(polynomial)};'
               'ideal J=std(jacob(F)+ideal('+','.join(['F', *map(str, earlier)])+'));'
               'print(dim(J));print(vdim(J));ideal low=std(J+ideal('+','.join(map(fmt, minors))+'));'
               'print(vdim(low));matrix H=jacob(jacob(F));print(vdim(std(J+ideal(det(H)))));')
    if interpolate:
        cubics = tuple(combinations_with_replacement(range(3), 3))
        for L, R in product(cubics, repeat=2):
            monomial = s.prod(l[h] for h in L)*s.prod(r[h] for h in R)
            program += f'print(string(reduce({fmt(monomial.subs(fixed))},J)));'
    result = subprocess.run([os.getenv('SINGULAR_BIN','Singular'), '-q'], input=program+'quit;',
                            text=True, capture_output=True, check=True)
    lines = [z for z in result.stdout.splitlines() if z and not z.startswith('//')]
    assert len(lines) == 4+100*interpolate and '?' not in result.stdout, result.stdout
    dim, size, low, degenerate = map(int, lines[:4])
    assert dim <= 0 and degenerate == 0
    values = [s.Poly(s.sympify(line.replace('^','**')), *variables).as_dict() for line in lines[4:]]
    monomials = sorted(set().union(*(set(f) for f in values)))
    return size, low, [[int(f.get(m,0))%PRIME for f in values] for m in monomials]


def nodal_counts(data):
    l, r = s.symbols('a b c'), s.symbols('d e f')
    G, F, Q, P = data
    for first, second, interpolate in ((Q,P,True), (G,F,False)):
        C = bracket(first, second)
        M = s.Matrix(3, 3, lambda a,b: sum(int(C[i,j,a,b])*l[i]*r[j] for i,j in product(range(3),repeat=2)))
        size = low = 0
        values = []
        for i,j in product(range(3), repeat=2):
            n, m, samples = checked_chart(M,l,r,i,j,interpolate)
            size += n
            low += m
            values.extend(samples)
        assert (size,low) == (71,36)
        if interpolate:
            assert nmod_mat(values, PRIME).rank() == 68


def tangent(data):
    p = PRIME
    def polarized(left, right):
        g, _, q, _ = left
        _, f, _, p0 = right
        return (bracket(g,f)+bracket(q,p0).transpose(2,3,0,1)).ravel()%p
    assert not any(polarized(data,data))
    slots = []
    for t, i, k in product(range(4),range(3),range(8)):
        a, b = divmod(k,3)
        slot = np.zeros_like(data)
        slot[t,i,a,b] = 1
        if a == b:
            slot[t,i,2,2] = -1
        slots.append(slot)
    J = np.array([(polarized(v,data)+polarized(data,v))%p for v in slots]).T
    basis = kernel(rows(J,p))
    assert len(basis) == 34
    directions = (basis@np.array(slots).reshape(96,-1)%p).reshape(-1,4,3,3,3)
    rng = np.random.default_rng(345)
    values = []
    for _ in range(130):
        l,r = rng.integers(p, size=(2,3))
        q = (l@data[2].reshape(3,9)%p).reshape(3,3)
        p0 = (r@data[3].reshape(3,9)%p).reshape(3,3)
        M = (q@p0-p0@q)%p
        cof = np.array([np.cross(M[(h+1)%3], M[(h+2)%3])%p for h in range(3)])
        def differential(dq,dp):
            dM = (dq@p0+q@dp-dp@q-p0@dq)%p
            return np.sum(cof*dM, axis=(-1,-2))%p
        zero = np.zeros((3,3), dtype=np.int64)
        dl = [differential(m,zero) for m in data[2]]
        dr = [differential(zero,m) for m in data[3]]
        orbit = [*np.outer(l,dl).ravel()%p, *np.outer(r,dr).ravel()%p, M[0]@cof[0]%p]
        motion = differential(np.einsum('j,hjbc->hbc',l,directions[:,2])%p,
                              np.einsum('j,hjbc->hbc',r,directions[:,3])%p)
        values.append([*orbit, *motion])
    V = np.array(values)
    assert [rows(m,p).rank() for m in (V[:,:19],V[:,19:],V)] == [17,24,32]


if __name__ == '__main__':
    triples = marked_triples()
    nodal_counts(triples)
    tangent(triples)
    print('PASS: 71 ordinary nodes in each direction, split 36+35, with evaluation rank 68 and differential rank 32')
