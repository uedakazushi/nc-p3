"""Parallel masks 94--96, their two surfaces and their unmarked parameter quotients."""
from itertools import combinations, product

import sympy as s

from verify_klein_mixed_positive import bit
from verify_klein_sparse_transverse import torus_forms
from verify_klein_three_lines import K, twist
from verify_klein_two_doublet_supports import candidates
from verify_transverse_rank_sixteen_isomorphism import adjacent, opposite, universal
from verify_transverse_rank_twelve_isomorphism import alternating_gradient, radical_determinant


def symmetric(p,q):
    values = {(0,0,0,1):p,(0,1,1,1):1,(0,1,2,2):1,(0,0,2,3):1,(2,2,2,3):q,(2,3,3,3):1}
    return {k:s.sympify(values.get(tuple(sorted(k)),0)) for k in product(K,repeat=4)
            if k[0]^k[1]^k[2]^k[3]==1}


def supports(p,q):
    for n,P,Q in ((94,0,0),(95,p,0),(96,p,q)):
        raw = {(i,j^1,k,l):v for (i,j,k,l),v in symmetric(P,Q).items() if v}
        mask = candidates(1)[n]
        assert sum(bit(t) for t in raw)==mask
        assert [(len(y),eq) for _,y,eq in torus_forms(mask,1)]==[(n-94,[])]*2
        gauge = s.Matrix([[int(t[h]==i) for h,i in product(K,repeat=2)] for t in raw])
        parameters = s.Matrix([[int(v==z) for z in (p,q)] for v in raw.values()])
        assert gauge.rank()==13 and gauge.row_join(parameters).rank()==13+n-94
    exponents = s.Matrix([[1,3,0,0],[1,1,2,0],[2,0,1,1],[0,0,1,3]])
    extra = s.Matrix([[3,1,0,0],[0,0,3,1]])
    assert 2*extra*exponents.inv()==s.Matrix([[1,-1,3,-1],[-1,3,-1,1]])


def web(S,p,q):
    # Coordinates X,Y,Z,T. Rows are the four quadratic equations for a square in R.
    normals = [[(1,-p,0,-1),(0,-1,1,-q)],[],[(0,1,0,1)],[(1,-1,-1,1)]]
    for h in K:
        B = adjacent(S,h)
        assert B.rank()==(2 if h<2 else 1)
        for row in normals[h]:
            assert s.Matrix([row])*B==s.zeros(1,4)
    a,b,u,t = s.symbols('a b u t')
    M = s.Matrix([[a,0,0,u],[0,-p*a-b,-u,t],[0,-u,b,0],[u,t,0,-a-q*b]])
    assert s.factor(M.det().subs(t,0)-((p*a+b)*b+u**2)*((a+q*b)*a+u**2))==0
    f = b**2+(p-q)*b-1
    assert s.discriminant(f,b)==(p-q)**2+4 or s.simplify(s.discriminant(f,b)-(p-q)**2-4)==0
    assert s.factor(s.discriminant(f,b).subs(q,1/p)-(p+1/p)**2)==0
    # The restriction to the kernel at a=0, b=1, u=0, t=t0, t0^2=q.
    t0 = s.Symbol('t0')
    H = s.Matrix([[1,0],[0,t0],[0,0],[0,1]])
    R = H.T*M.subs(q,t0**2)*H
    assert s.expand(R[1,1]+(1+p*t0**2)*a+2*t0**2*b-2*t0*t)==0
    assert R[0,0]==a and R[0,1]==u


def covariant(S,p,q):
    g = alternating_gradient(S)
    D = s.factor(sum(S[t]*g[t] for t in S)/4)
    z = p*q
    assert D==8*(2*z+3) or s.expand(D-8*(2*z+3))==0
    L = [(adjacent(S,h)*adjacent(g,h).T).applyfunc(s.factor) for h in K]
    v = [s.factor(sum((B**3)[i,i] for B in L)) for i in K]
    assert v[0]==v[2] and v[1]==v[3] and s.factor(v[0]-v[1]-2048*z*(z-1))==0
    T = s.factor(sum((B**2).trace() for B in L)/D**2)
    assert T.subs(p,0)==s.Rational(8,3) and s.cancel(T.subs(q,1/p))==s.Rational(88,25)


def twisted(S,p,q):
    W = {t:v*twist(t) for t,v in S.items()}
    for h in (0,1):
        assert opposite(W,h).det()==-4
    expected = [[(0,1,0,0),(0,0,0,1)],[(1,0,1,0),(0,-1,0,1)]]
    for h in (2,3):
        B = opposite(W,h)
        ker = [tuple(v) for v in B.T.nullspace()]
        assert ker==expected[h-2] and [tuple(v) for v in B.nullspace()]==expected[3-h]
    z = s.symbols('z:4')
    assert radical_determinant(W)==-z[2]**2*z[3]**2
    boundary = universal(1,p,0,1)
    assert all(S[t].subs(q,p)==boundary[tuple(i^1 for i in t)]*twist(tuple(i^1 for i in t)) for t in S)


if __name__=='__main__':
    p,q = s.symbols('p q')
    S = symmetric(p,q)
    supports(p,q)
    web(S,p,q)
    covariant(S,p,q)
    twisted(S,p,q)
    print('PASS: six torus components, all square equations, intrinsic splitting, rank-two marks and transverse boundary')
