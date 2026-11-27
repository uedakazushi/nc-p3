"""Fourier equations and regular strata at the simultaneous full-thin degeneration."""
from itertools import combinations, product

import sympy as s

from verify_negative_component_quotient import K, thin, words

g0,g1,r0,r1 = s.symbols('g0 g1 r0 r1', nonzero=True)
U,V = r0+r1,g0+g1
kappa = 1-U*V


def twist(q):
    return (-1)**sum((q[i]>>1)*(q[j]&1) for i,j in combinations(K,2))


def coefficients():
    extra = {(0,0,0,0):kappa/g0,(1,1,1,1):kappa/g1,
             (2,3,2,3):kappa/r0,(3,2,3,2):kappa/r1}
    return {q: (1 if q in thin else (U if max(q)<2 else V)+extra.get(q,0)) for q in words}


def check_thin_twist():
    # Four one-letter signs turn the cocycle exactly into lambda -> -lambda on thin words.
    factors = [[-1,-1,1,1],[-1,1,1,1],[1,-1,-1,1],[1,1,-1,1]]
    for q in words:
        i,j,k,_ = q
        h = i^j
        sign = twist(q)*s.prod(factors[p][v] for p,v in enumerate(q))
        if q in thin:
            expected = (-1 if h==0 and (i^k)&1 or h==1 and not (i^k)&1 else 1)
            assert sign==expected
        elif max(q)<2:
            assert sign==(-1)**(h+(i&1)+(k&1))


def fourier():
    m,n,A,B,C,D,E,F = s.symbols('m n A B C D E F')
    matrix = s.Matrix([[m,B,1,0],[A,C,0,0],[1,0,n,E],[0,0,D,F]])
    minors = {s.factor(matrix.extract(i,j).det()) for i in combinations(K,3) for j in combinations(K,3)}
    for v in (C*F,C*D,C*E,(m*C-A*B)*n-C):
        if v==(m*C-A*B)*n-C:
            assert any(s.expand(p-v)==0 or s.expand(p+v)==0 for p in minors)
        else:
            assert v in minors or -v in minors
    no_pairs = {p.subs({C:0,F:0}).expand() for p in minors} - {0}
    gens = {A*E,B*D,n*A*B,m*D*E}
    monomials = {p.as_coeff_Mul()[1] for p in no_pairs}
    assert gens <= monomials and all(any(s.denom(p/g)==1 for g in gens) for p in monomials)


def ranks_and_relations():
    W = coefficients()
    def opposite(w,h):
        return s.Matrix(4,4,lambda i,j:w[i,j,i^h,j^h])
    first = [s.Matrix(4,4,lambda i,j:W[i,i^h,j,j^h]) for h in K]
    targets = [s.Matrix([[1,0],[0,1],[g0,g1],[g0,g1]]),
               s.Matrix([[r0,r1],[r0,r1],[1,0],[0,1]])]
    for h,target in enumerate(targets):
        select = [0,1] if h==0 else [2,3]
        assert (first[h]-target*first[h][select,:]).applyfunc(s.factor)==s.zeros(4)
        assert s.factor(first[h].extract(select,select).det())!=0
    assert all(first[h]==s.ones(4) for h in (2,3))
    assert s.factor(opposite(W,0).det()*g0*g1*r0*r1/kappa**3)==-1
    assert opposite(W,1)[:2,2:]==s.ones(2) and opposite(W,1)[2:,:2]==s.ones(2)
    for h in (2,3):
        assert opposite(W,h)==s.ones(4)
    W = {q:v*twist(q) for q,v in W.items()}
    disc = kappa**2+(4*U*g0*g1+kappa*V)*(4*V*r0*r1+kappa*U)
    for h,det in enumerate((-disc/(g0*g1*r0*r1),-16*U*V,-16,-16)):
        assert s.factor(opposite(W,h).det()-det)==0
    # Solvability at an all-nonzero point would force U*V=1.
    a,b=s.symbols('a b',nonzero=True)
    elimination = (r1/b-r0*a)*(g1/a-g0*b)-(-r0+r1*a/b)*(-g0+g1*b/a)
    assert s.factor(elimination-(g1/a+g0*a)*(r1/b+r0*b)+U*V)==0


def smoothing():
    a,b,d,j,e=s.symbols('a b d j e')
    l=1+e
    f=b+e*j
    t=(a+d-(l+1/l)*b)/(a*d-b*b)
    g=s.cancel((f+t*b*b-(l+1/l)*b)/(t*f-1))
    B=s.Matrix([[b,f],[g,b]])
    T=s.Matrix([[1,l],[l,1]])*B.inv()*s.Matrix([[1,1/l],[1/l,1]])
    v=(a+d-2*b)/(a*d-b*b)
    def residue(expr):
        return s.factor(s.cancel(expr).subs(e,0))
    assert residue(g)==b and residue(T[0,1]-v)==residue(T[1,0]-v)==0
    for i,den in enumerate((b+j,b-j)):
        assert s.factor(residue(T[i,i]-v)*den-2*(1-b*v))==0
    assert residue(B.det()/e**2)!=0


if __name__ == '__main__':
    check_thin_twist()
    fourier()
    ranks_and_relations()
    smoothing()
    print('PASS: thin sign change, all Fourier pivot equations, normalized relations, ranks, point elimination and smoothing')
