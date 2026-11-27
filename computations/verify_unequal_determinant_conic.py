"""Check the small-minor obstruction, inverse point fibers and rational quotient."""
from itertools import combinations, product

import sympy as S

from verify_exchange_mixed_configurations import incidence_signatures, maps, normalize, partner
from verify_exchange_singular_corners import index_row
from verify_full_thin_exchange_cover import tensor
from verify_klein_elliptic_completion import gram_twist
from verify_negative_component_quotient import K, relabel, thin, words
from verify_parallel_invertible_rectangles import smith

a,b,d,k = S.symbols('a b d k')
ratio = (k*k-b*b)/((a-k)*(d-k))
D = a*d-b*b
e,h = k+ratio*(a-k),k+ratio*(d-k)
u,v,s,t = x = S.symbols('u v s t')
m,n,f,g,M,N,F,G = coeff = S.symbols('m n f g M N F G')
sub = dict(zip(coeff,[(d-b)/D,(d+b)/D,(a-b)/D,(a+b)/D,
                      (h-b)/(ratio*D),(h+b)/(ratio*D),(e-b)/(ratio*D),(e+b)/(ratio*D)]))
E = S.Matrix([[u,0,m*s,-n*t],[0,v,f*s,g*t],[0,u,N*t,M*s],[v,0,-G*t,F*s],
              [s,-t,u,v],[t,s,v,u]])


def identities():
    minor = lambda I:E.extract(I,K).det(method='domain-ge')
    poly = S.Poly((m-M)*minor((0,1,2,3))+M*f*minor((0,3,4,5))-F*m*minor((1,2,4,5)),*x)
    for exponent,value in poly.terms():
        value = S.factor(value.subs(sub))
        if exponent==(1,1,1,1):
            expected = -2*(a-b)*(a+b)*(b-d)*(b+d)*(a-k)*(d-k)*(1+ratio)/((b-k)*(b+k)**2*D**3)
            assert S.factor(value-expected)==0
        else:
            assert value==0
    A = S.Matrix([[a,b],[b,d]])
    B = S.Matrix([[e,b],[b,h]])
    w,op,_ = gram_twist(A,B,-S.S.One,S.S.One)
    rows = [S.Matrix(4,4,lambda i,j:S.expand(E[r,j]).coeff(x[i])).subs(sub).reshape(16,1) for r in range(6)]
    # Compare the column spaces in each color block.
    for j,indices in enumerate(((0,1),(2,3))):
        ix = [4*i+(i^j) for i in K]
        U = S.Matrix.hstack(*(rows[r] for r in indices))[ix,:]
        assert U[:2,:].det() in (1,-1)
        assert (w[ix,ix]-U*U[:2,:].inv()*w[ix[:2],ix]).applyfunc(S.factor)==S.zeros(4)
    for j in (2,3):
        ix = [4*i+(i^j) for i in K]
        col = rows[j+2][ix,:]
        assert (w[ix,ix]-col*w[ix[:1],ix]).applyfunc(S.factor)==S.zeros(4)
    block = op[[0,5,10,15],[0,5,10,15]].applyfunc(S.factor)
    det = S.factor(block.det(method='domain-ge'))
    factors = S.factor_list(S.fraction(det)[0])[1]
    expected = [(S.factor((a-k)*(d-k)*(1-ratio)),2),(S.factor((a-k)*(d-k)*(1+ratio)),1)]
    assert set(factors)==set(expected)
    for j in (1,2,3):
        ix = [4*i+(i^j) for i in K]
        block = op[ix,ix].applyfunc(S.factor)
        assert block.rank()==2
        if j==1:
            assert block[0,2:4]==S.Matrix([[1,1]]) and block[2,:2]==S.Matrix([[1,-1]])
            assert block[0,0]==block[0,1]
        else:
            assert any(S.factor(block.extract(I,J).det()) in (2,-2)
                       for I in combinations(K,2) for J in combinations(K,2))
    Z2 = -m*N/(n*M)
    X2 = m+f*m/M
    # y=(X,M*Z*X/m,1,Z) has rank two on ratio=-1.
    for expr in (F*N+f*n,G*M+g*m,n-g*m/M-M*X2/m,
                 M*X2*Z2/m-m*Z2-f*m/M,X2+n*Z2-g*m/M):
        assert S.factor(S.factor(expr.subs(sub)).subs(d,k+(b*b-k*k)/(a-k)))==0


def image_counts():
    A = S.Matrix([[5,2],[2,11]])
    B = S.Matrix([[7,2],[2,37]])
    w,_,_ = gram_twist(A,B,-S.S.One,S.S.One)
    p = 761
    rows = [(1,1,0,0),(0,0,1,-1),(0,-96,-507,315),(-24,0,-63,15),
            (-507,-63,0,3774),(315,15,2550,0)]
    pts = set()
    for q in rows:
        c = next(i for i in q if i)
        roots = [S.sqrt_mod(i*pow(c,-1,p)%p,p,all_roots=True) for i in q]
        pts.update(normalize(z,p) for z in product(*roots))
    assert len(pts)==20
    for config,order in ((sorted(pts),8),(partner(sorted(pts),w,False,p),4)):
        hist,sig = incidence_signatures(config,p)
        assert max(hist)==6
        assert maps((config,sig),(config,sig),p)==order


def marked_action():
    a,b,d,e,h = S.symbols('a b d e h')
    A,B = S.Matrix([[a,b],[b,d]]),S.Matrix([[e,b],[b,h]])
    T = tensor(A,B,-S.S.One,S.S.One)
    top = [(0,0,0,0),(0,0,1,1),(1,1,1,1),(0,1,0,1),(1,0,1,0)]
    W = S.Matrix([index_row(q) for q in thin])
    H,L,R = smith(W)
    rank = W.rank()
    found = set()
    for shear,i in product(range(2),K):
        for j,k in product((i&2,(i&2)^1),repeat=2):
            args = 0,shear,(i,j,k,i^j^k)
            V = {q:T[relabel(q,*args)] for q in words}
            rhs = L*S.Matrix([2*int(V[q]!=T[q]) for q in thin])
            if any(v%4 for v in rhs[rank:]):
                continue
            phase = R*S.Matrix([rhs[t]/H[t,t] if t<rank else 0 for t in range(16)])
            vec = [V[q]*S.I**(int(sum(phase[4*t+m] for t,m in enumerate(q)))%4) for q in top]
            aa,bb,dd,ee,hh = vec
            found.add(tuple(S.factor(z/bb) for z in (aa,dd,ee,hh)))
    beta = (e*h-b*b)/(a*d-b*b)
    expected = set()
    for aa,dd,ee,hh in [(a/b,d/b,e/b,h/b),
                         ((2*b-a-d)/(a-d),(-2*b-a-d)/(a-d),
                          (2*b-e-h)/(beta*(a-d)),(-2*b-e-h)/(beta*(a-d))),
                         ((2*b-e-h)/(e-h),(-2*b-e-h)/(e-h),
                          beta*(2*b-a-d)/(e-h),beta*(-2*b-a-d)/(e-h))]:
        for z in ((aa,dd,ee,hh),(dd,aa,hh,ee),(ee,hh,aa,dd),(hh,ee,dd,aa)):
            for sign in (1,-1):
                expected.add(tuple(S.factor(sign*v) for v in z))
    # The conic identity e-h=beta*(a-d) identifies the two denominators.
    hh = (e*a*d-e*b*b+(a-d)*b*b)/(a*d-b*b+(a-d)*e)
    restrict = lambda values:{tuple(S.factor(z.subs(h,hh)) for z in v) for v in values}
    assert restrict(found)==restrict(expected) and len(restrict(found))==16


if __name__=='__main__':
    identities()
    image_counts()
    marked_action()
    print('PASS: quartic identity, bad ratio, opposite ranks, two frame stabilizers and all marked actions')
