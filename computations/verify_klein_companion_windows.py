"""Check the Klein flatness equations and an independent nontransverse PBW example."""
from itertools import combinations, product

import sympy as sp
from verify_twenty_point_companion import I, pairs


def verify_simple_equations():
    alpha, delta, b, c = [sp.symbols(f'{name}:3') for name in ('a','d','b','c')]
    a, bs, cs = [], [], []
    for i in range(3):
        j,k = (i+1)%3, (i+2)%3
        m = sp.zeros(3)
        m[k,j],m[j,k] = alpha[i],delta[i]
        a.append(m)
        t = sp.zeros(3)
        t[j,k] = t[k,j] = 1
        bs.append(b[i]*t)
        cs.append(c[i]*t)
    def value(x,y,dual):
        u,v = x[1:,:],y[1:,:]
        ops = [m.T for m in a] if dual else a
        left = sum((u[i]*ops[i] for i in range(3)), sp.zeros(3))
        right = sum((v[i]*m for i,m in enumerate(cs if dual else bs)), sp.zeros(3))
        return (x[0]*v+left*v).col_join(y[0]*u+right*u)
    mats = [sp.Matrix.hstack(*(value(I[:,i],I[:,j],d) for i,j in pairs)) for d in (False,True)]
    w = mats[0].T*mats[1]
    errors = {sp.expand(w[4*i+j,4*k+l]-w[4*k+j,4*i+l])
              for i,j,k,l in product(range(4),repeat=4)}
    eq = {b[i]*c[(i+1)%3]+delta[i]*delta[(i+1)%3] for i in range(3)}
    eq |= {c[i]*b[(i+1)%3]+alpha[i]*alpha[(i+1)%3] for i in range(3)}
    assert errors == {0} | eq | {-v for v in eq}


def matching_tensor(m):
    square = sp.zeros(4)
    for h in (1,2,3):
        a,b = 0,h
        c,d = [i for i in range(4) if i not in (0,h)]
        p = m[[b,a],[d,c]]
        if h == 1:
            assert p.rank() == 1
            continue
        x,y,z,t = list(p)
        u = sp.diag(1,-z*t/(x*y))
        v = -p.T*u.inv()*p
        assert v[0,1] == v[1,0] == 0
        square[a,b],square[b,a] = u[0,0],u[1,1]
        square[c,d],square[d,c] = v[0,0],v[1,1]
    def coefficient(a,b,c,d):
        if a == c and b == d:
            return square[a,b]
        return m[b,d] if len({a,b,c,d}) == 4 else 0
    w = sp.Matrix([[coefficient(a,b,c,d) for c,d in pairs] for a,b in pairs])
    return w,square


if __name__ == '__main__':
    verify_simple_equations()
    m = sp.zeros(4)
    for (a,b),v in zip(combinations(range(4),2),[2,1,1,3,3,5]):
        m[a,b],m[b,a] = v,-v
    w,d = matching_tensor(m)
    assert d.det() == 36
    opp = sp.Matrix([[w[4*a+b,4*c+e] for b,e in pairs] for a,c in pairs])
    assert w.rank() == 6 and opp.rank() == 10
    for h in (1,2,3):
        ids = [4*a+(a^h) for a in range(4)]
        assert w[ids,ids].rank() == 2
    descents = [4*i+j for i,j in pairs if i>j]
    for a in [w,w.T]:
        rel = sp.Matrix.hstack(*a.columnspace()).T
        assert rel[:,descents].rank() == 6
        assert rel[:,[4*(j%4)+j//4 for j in descents]].rank() == 6
    print("Six symbolic flatness equations; matched PBW example, adjacent rank 6 and opposite rank 10.")
