"""Exclude every simultaneous (3,1,1,1) and (2,2,1,1) cut."""
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path
from subprocess import check_output, run
from tempfile import TemporaryDirectory

import sympy as s
from sympy.polys.matrices import DomainMatrix

from verify_klein_three_lines import K, WORDS, character, entry, smith, thin_system
from verify_klein_single_missing_supports import admissible, upper

MASKS = (0xff00f0f00ffcfff,0xfff0fff0fffcfff)


def bit(w):
    i,j,k,_ = w
    return 1<<(16*(i^j)+4*i+k)


def rank(p,b):
    return (1 if b else 3) if p%2==0 else (2 if b<2 else 1)


def sparse_supports():
    source = Path(__file__).with_name('klein_mixed_positive_supports.cpp')
    with TemporaryDirectory() as temp:
        binary = Path(temp)/'supports'
        run(['c++','-std=c++17','-O3',str(source),'-o',str(binary)],check=True)
        masks = set(map(int,check_output([str(binary)],text=True).splitlines()))
    assert len(masks)==14744
    @lru_cache(None)
    def bound(mask,r,indices):
        if r<3:
            return upper(mask,r,indices)
        return min(r,max(sum(mask>>(4*i+pi[i])&1 for i in indices) for pi in permutations(K)))
    n = 0
    surviving = []
    while masks:
        mask = next(iter(masks))
        words = [w for w in WORDS if mask&bit(w)]
        orbit = set()
        for seq in (words,[w[2:]+w[:2] for w in words],[(b,a,d,c) for a,b,c,d in words],
                    [w[::-1] for w in words]):
            for pi in ((0,1,2,3),(0,1,3,2)):
                for t,h in product(K,(0,1)):
                    orbit.add(sum(bit((pi[a]^t,pi[b]^t,pi[c]^t^h,pi[d]^t^h)) for a,b,c,d in seq))
        assert orbit<=masks
        masks -= orbit
        n += 1
        if admissible(words,rank,lambda m,r,I:bound(m,r,tuple(I))):
            surviving.append(min(orbit))
    assert n==312 and tuple(sorted(surviving))==MASKS


def normal_form(words):
    thin,ds,right,gauge = thin_system([(1,2,3),(2,3)],words)
    assert ds==[1]*len(ds)
    _,left,_ = smith(gauge)
    n = gauge.rows-gauge.rank()
    coords = [s.Integer(1)]*gauge.rank()+list(s.symbols('a:'+str(n)))
    free = [character(row,coords) for row in left.inv().tolist()]
    base = [s.Integer(1)]*len(ds)+free
    w = {key:character(row,base) for key,row in zip(thin,right.tolist())}
    missing = sorted(set(words)-set(thin))
    assert len(missing)==8
    w.update(zip(missing,s.symbols('x:8')))
    return w


def matrix(w,p,b):
    return s.Matrix([[w.get(entry(p,b,i,k),0) for k in K] for i in K])


def minors(m):
    return [s.fraction(s.factor(m[I,J].det()))[0] for I,J in product(combinations(K,3),repeat=2)]


def dense():
    w = normal_form(WORDS)
    a = s.Symbol('a0')
    x = s.symbols('x:8')
    z = s.Symbol('z')
    for b,indices in [(0,(0,3,4,7)),(1,(1,2,5,6))]:
        m = matrix(w,1,b)
        var = [x[i] for i in indices]
        u,v,f,h = var
        def localized(polys):
            return list(s.groebner(polys+[z*a*(a*a-1)-1],*var,z,a,order='grevlex',domain=s.QQ))
        eq = minors(m)
        assert localized(eq)==localized([u-v,u*h-1,f-h,(h-a)*(a*h-1)])
        for val in (1,-1):
            left = [q.subs(a,val) for q in eq]
            pairs = [(j-val)*(k-val) for j,k in combinations(var,2)]
            assert list(s.groebner(left,*var))==list(s.groebner(pairs,*var))
            assert m.subs({a:val,**dict.fromkeys(var,val)}).rank()==1
    m = matrix(w,0,0)
    for b,c in product((a,1/a),repeat=2):
        out = m.subs(dict(zip(x,(1/b,1/c,1/c,1/b,b,c,c,b))))
        if b==c:
            assert out.rank()==1
        else:
            assert s.factor(out.det()+(a-1/a)**4)==0
    for val in (1,-1):
        assert m.subs(dict.fromkeys(x,val)).rank()==1


def fourth_character(w,p,g):
    cols = [(i,j,k,i^j^k^g) for i,j,k in product(K,repeat=3)]
    rows = []
    for t,b in product(range(3),K):
        mat = matrix(w,(p+t)%4,b)
        for col in [mat[:,k] for k in K]:
            for j in K:
                outer = (j,j^b^g)
                v = {outer[:t]+(i,i^b)+outer[t:]:val for i,val in enumerate(col) if val}
                rows.append([v.get(word,0) for word in cols])
    return 64-DomainMatrix.from_Matrix(s.Matrix(rows)).to_field().rank()


def sparse_equations(n,mask):
    w = normal_form([key for key in WORDS if mask&bit(key)])
    a,b = s.symbols('a:2')
    x = s.symbols('x:8')
    eq = minors(matrix(w,1,0))+minors(matrix(w,1,1))
    # Work in the nonzero coefficient torus. Remove only monomial unit factors.
    def primitive(f):
        poly = s.Poly(f,a,b,*x,domain=s.QQ)
        return poly.terms_gcd()[1].monic().as_expr() if f else 0
    equations = set(map(primitive,eq))
    required = ([a-b,x[0]-b,a*x[3]-1,b*x[4]-a,a**3*x[0]-1,
                 x[1]-1,x[2]-1,x[5]-1] if n==0 else
                [x[0]-a,x[1]-a,x[2]-1,x[3]-1,x[4]-a,x[5]-a,a*x[1]-1])
    if n==0:
        # The two equations x1=1 and a*x1=b together imply a=b.
        required[0] = a*x[1]-b
    assert set(map(primitive,required))<=equations
    vals = (1,-1,s.I,-s.I) if n==0 else (1,-1)
    for value in vals:
        t = s.sympify(value)
        sub = {a:t,b:t}
        first = (t,1,1,1/t,1,1) if n==0 else (t,t,1,1,t,t)
        sub.update(zip(x[:6],map(s.sympify,first)))
        if t in (1,-1):
            assert matrix(w,0,0).subs(sub).rank()<=2
        else:
            assert s.expand(matrix(w,0,0).subs(sub).det()+4*(x[6]+x[7]))==0
            sub[x[7]]=-x[6]
            candidate = {key:v.subs(sub) for key,v in w.items()}
            assert fourth_character(candidate,1,2)==9


if __name__=='__main__':
    sparse_supports()
    dense()
    for n,mask in enumerate(MASKS):
        sparse_equations(n,mask)
    print('PASS: all mixed positive profiles excluded, including 312 sparse orbits and all exceptional values')
