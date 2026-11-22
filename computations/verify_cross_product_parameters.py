"""Polynomial identities for the diagonal three-axis family (all ranks over characteristic zero)."""
from collections import Counter
from itertools import combinations, product

import sympy as sp
from verify_twenty_point_companion import I, J, flip, pairs, tensor


def marked_matrix(parameters):
    scaling = sp.diag(*parameters)
    def value(x, f):
        u, v = x[1:, :], f[1:, :]
        t = scaling*v
        mixed = sp.Matrix([u[(i+1)%3]*t[(i+2)%3]+u[(i+2)%3]*t[(i+1)%3]
                           for i in range(3)])
        return (x[0]*v+(scaling*u).cross(v)).col_join(f[0]*u+mixed)
    return sp.Matrix.hstack(*(value(I[:,i], I[:,j]) for i,j in pairs))


def parameters(lam, parity):
    rows = [[0,*lam]]
    for i in range(3):
        j, k = (i+1)%3, (i+2)%3
        row = [-lam[i],0,0,0]
        row[j+1] = (-1)**parity*lam[k]*(lam[i]**2+lam[j]**2)
        row[k+1] = -(-1)**parity*lam[j]*(lam[i]**2+lam[k]**2)
        rows.append(row)
    return sp.Matrix(rows)


def verify():
    lam = sp.symbols('a b c')
    a,b,c = lam
    delta = a*a*b*b+a*a*c*c+b*b*c*c
    matrix = marked_matrix(lam)
    companion = -matrix*tensor(J,J)
    w = matrix.T*companion
    assert all(sp.expand(w[4*i+j,4*k+l]-w[4*k+j,4*i+l]) == 0
               for i,j,k,l in product(range(4),repeat=4))
    def block(g):
        rows = [(i,k) for i,k in pairs if i <= k and i^k == g]
        cols = [(j,l) for j,l in pairs if j^l == g]
        return sp.Matrix([[w[4*i+j,4*k+l] for j,l in cols] for i,k in rows])
    assert sp.factor(block(0).det()+4*delta) == 0
    for i in range(3):
        m = block(i+1)
        assert m.rank() == 2 and all(sp.cancel(t/lam[i]).is_Integer for t in m)
        assert any(sp.expand(m[:,ks].det()) in [lam[i]**2,-lam[i]**2]
                   for ks in combinations(range(4),2))
    rels = [matrix, matrix[:,flip]*tensor(I,J)]
    for p in range(2):
        q, q_next = parameters(lam,p), parameters(lam,1-p)
        assert sp.factor(q.det()) == 4*delta**2
        relations = sp.Matrix.vstack(tensor(rels[p],I),tensor(I,rels[1-p]))
        for i in range(4):
            words = [(j,k,i^j^k) for j,k in pairs]
            rr = relations[:,[16*j+4*k+l for j,k,l in words]]
            rr = rr[[j for j in range(rr.rows) if any(rr[j,:])],:]
            identities = []
            for h in range(4):
                chi = 1 if h==0 or i in (0,h) else -1
                row = [0]*len(words)
                for j in range(4):
                    row[words.index((j,j,i))] += q[h,j]
                    row[words.index((i,j,j))] -= chi*q_next[h,j]
                identities.append(row)
            assert rr.to_DM().rank() == rr.col_join(sp.Matrix(identities)).to_DM().rank() == 11
    print("Verified flatness, opposite blocks, Pfaffian and all signed normalities over Q(a,b,c).")


def check_trisecants():
    # A reduced split configuration over F_97, checked directly from all triples.
    prime = 97
    def root(t):
        return next(x for x in range(prime) if (x*x-t) % prime == 0)
    def normalize(v):
        h = pow(next(x % prime for x in v if x % prime), -1, prime)
        return tuple(h*x % prime for x in v)
    points = [tuple(int(i==j) for j in range(4)) for i in range(4)]
    h = root(-pow(2,-1,prime))
    for i in range(3):
        for u,v in product([1,-1], repeat=2):
            tail = [0]*3
            tail[(i+1)%3], tail[(i+2)%3] = u*h, v*h
            points.append(normalize([1,*tail]))
    points.extend(normalize([0,1,u,v]) for u,v in product([1,-1],repeat=2))
    assert len(set(points)) == 20
    counts = Counter()
    for i,j in combinations(range(20),2):
        v,w = points[i],points[j]
        line = {normalize([a+t*b for a,b in zip(v,w)]) for t in range(prime)} | {w}
        ids = [k for k,z in enumerate(points) if z in line]
        if len(ids)>2 and ids[:2]==[i,j]:
            assert len(ids)==3
            counts.update(ids)
    assert [counts[i] for i in range(20)] == [6]*4+[3]*16


if __name__ == '__main__':
    verify()
    check_trisecants()
