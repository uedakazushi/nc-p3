"""Rank-two relation curves, the marked binary graph and twisted double-pair separation."""
from itertools import combinations, product

import sympy as s

from verify_klein_single_missing_windows import regular_window
from verify_transverse_rank_sixteen_isomorphism import adjacent, universal
from verify_klein_three_lines import K, WORDS, twist
from verify_transverse_rank_twelve_isomorphism import alternating_gradient


def relation_curves():
    x,y,z = s.symbols('x y z')
    diagonals = s.symbols('d:4')
    M = s.Matrix([[diagonals[i] if i==j else (-1)**((i>>1)*((i^j)&1))*(x,y,z)[(i^j)-1]
                   for j in K] for i in K])
    adj = M.adjugate()
    f = x**2+y**2-z**2
    for i,j in combinations(K,2):
        b = i^j
        k,l = [h for h in K if h not in (i,j)]
        sign = -1 if b==3 else 1
        symmetric = s.factor((adj[i,j]+sign*adj[j,i])/2)
        target = (x,y,z)[b-1]*(f-diagonals[k]*diagonals[l])
        assert s.factor(symmetric**2-target**2)==0
        skew = s.factor((adj[i,j]-sign*adj[j,i])/2)
        target = s.prod(t for n,t in enumerate((x,y,z),1) if n!=b)*(diagonals[k]-diagonals[l])
        assert s.factor(skew**2-target**2)==0


def exceptional():
    r,X = s.symbols('r X')
    a,b = (1+r)/4,(1-r)/4
    d0,d1 = s.symbols('d0 d1')
    w = {t:twist(t)*(1+((d0,d1)[t[0]//2] if len(set(t))==1 else 0)) for t in WORDS}
    g = alternating_gradient(w)
    W,G = [{t:s.factor(f.subs({d0:-1/a,d1:-1/b})) for t,f in v.items()} for v in (w,g)]
    cap = s.factor(sum(W[t]*G[t] for t in WORDS)/4)
    assert s.factor(cap-64*(3*r**4-5*r**2+6)/(r**2-1)**2)==0
    def block(tensor,b):
        return s.Matrix([[tensor[i,i^b,k,k^b] for k in K] for i in K])
    L = [(block(W,h)*block(G,h).T/cap).applyfunc(s.factor) for h in K]
    diffs = []
    for n in (2,3):
        F = s.factor(sum((B**n).trace() for B in L))
        diffs.append(s.fraction(s.factor(F-F.subs(r,1/r)))[0])
    assert s.factor(s.gcd(*diffs)/(r**4-1)).is_Rational
    assert s.gcd(3*X**2-5*X+6,6*X**2-5*X+3)==1
    A = s.Matrix([s.Symbol('u'),s.Symbol('v')])
    B = s.Matrix([[-a/b-1,1-a/b],[1-a/b,-a/b-1]])*A/2
    assert s.factor(sum(B)+(a/b)*sum(A))==0 and s.factor(B[0]-B[1]+A[0]-A[1])==0
    assert s.factor((a/b+1)**2/(1-a/b)**2-1/r**2)==0


def transverse_surface():
    p,q,u,v,y,z = s.symbols('p q u v y z')
    w = universal(p,q,1,1)
    sign = s.diag(1,1,1,-1)
    vectors = [((p,q,1,1),(1,1,p,q)),((1,1,0,0),(0,0,1,1)),
               ((1,1,1,1),),((1,1,-1,-1),)]
    for h,cols in enumerate(vectors):
        target = s.Matrix.hstack(*(s.Matrix(c) for c in cols))
        actual = sign*adjacent(w,h)
        assert target.row_join(actual).rank()==target.rank()==actual.rank()
    M = s.Matrix([[0,u,y,z],[u,0,z,y],[y,-z,0,v],[-z,y,v,0]])
    for f in M.adjugate():
        assert s.rem(f,u*v-y**2+z**2,u)==0


def missing_fixed_points():
    t = s.symbols('t')
    w0 = regular_window(t)
    windows = [{k:w0.get(k,s.Integer(0)) for k in WORDS},
               {k:s.Integer(twist(k))*(1+((2 if k[0]==0 else -2) if len(set(k))==1 else 0)) for k in WORDS}]
    def opposite(w,h):
        return s.Matrix([[w[i,j,i^h,j^h] for j in K] for i in K])
    z = s.symbols('z:4')
    ker = s.zeros(4)
    basis = iter(z)
    for h in K:
        B = opposite(windows[0],h)
        expected = (t**4+6*t**2+1 if h<2 else -(t**2-1)**2)*t**2/(t**2+1)
        assert s.factor(B.det()-expected)==0
        for c in B.subs(t,1).T.nullspace():
            var = next(basis)
            for i in K:
                ker[i,i^h]+=var*c[i]
    assert s.factor(ker.det()-(z[0]**2-z[2]**2)*(z[1]**2+z[3]**2)/2)==0
    for index,w in enumerate(windows):
        G = alternating_gradient(w)
        D = s.factor(sum(w[k]*G[k] for k in WORDS)/4)
        def adjacent_zero(W,h):
            return s.Matrix([[W[i,i^h,k,k^h] for k in K] for i in K])
        L = [(adjacent_zero(w,h)*adjacent_zero(G,h).T/D).applyfunc(s.factor) for h in K]
        T = s.factor(sum((m*m).trace() for m in L))
        if index==0:
            f = t**4+6*t**2+1
            assert s.gcd(s.prod(s.fraction(D)),f)==1
            num,den = s.fraction(s.factor(T-s.Rational(43,16)))
            assert s.rem(num,f,t)==0 and s.gcd(den,f)==1
        else:
            assert D==144 and T==s.Rational(124,27)


if __name__=='__main__':
    relation_curves()
    exceptional()
    transverse_surface()
    missing_fixed_points()
    print('PASS: all paired cofactors, binary graph and quartic-gradient separation on the exceptional line, transverse surface and singly missing intersections')
