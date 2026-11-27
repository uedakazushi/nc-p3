"""Simultaneous rank-two completions on both parallel thin shapes with invertible corners."""
from functools import lru_cache
from itertools import product

from flint import nmod_mat
import sympy as s

from verify_klein_mixed_positive import bit, matrix
from verify_klein_three_lines import K, WORDS, character, smith, thin_system
from verify_klein_two_doublet_supports import candidates
from verify_transverse_rank_twelve_isomorphism import alternating_gradient

T = s.Symbol('t',nonzero=True)
BETA,GAMMA,ALPHA = s.symbols('b c a',nonzero=True)
X = s.symbols('x:16')
A,B,C,D = s.symbols('A B C D')


@lru_cache(None)
def template(shape):
    mask = candidates(1)[shape]
    keys = [t for t in WORDS if (t[0]^t[1]<2 and t[1]^t[2]<2) or mask&bit(t)]
    thin,ds,R,H = thin_system([(2,3)]*2,keys)
    assert all(d==1 for d in ds)
    _,P,_ = smith(H)
    assert H.rows-H.rank()==1
    free = [character(row,[s.S.One]*H.rank()+[T]) for row in P.inv().tolist()]
    values = [character(row,[s.S.One]*len(ds)+free) for row in R.tolist()]
    w = dict(zip(thin,values))
    holes = [k for k in keys if k not in w]
    assert len(holes)==16
    w.update(zip(holes,X))
    return {k:s.sympify(v) for k,v in w.items()}


def top(shape,t,b,c,A,B,C,D):
    x0 = (c*A+(t**2-c)*B)/t if shape==76 else (c*A+t*(t**2-c)*B)/t
    x7 = c*C/t if shape==76 else ((t**2-c)*B+c*C)/t
    x3 = (t**2*D+(b-1)*c*C)/b
    if shape==112:
        x3 = (t**2*D+(b-1)*c*A+t*(b-1)*c*C+t*(t**2*(b-2)-2*b*c+b+2*c)*B)/b
    return tuple(map(s.factor,(x0,b*B,A,x3,B,C,D,x7)))


def completed(shape,t,b,c,A,B,C,D):
    w = template(shape)
    t,b,c,A,B,C,D = map(s.sympify,(t,b,c,A,B,C,D))
    sub = dict(zip(X,top(shape,t,b,c,A,B,C,D)))|{T:t}
    for h in (0,1):
        M = matrix(w,0,h).subs(sub)
        bottom = M[2:,:2]*M[:2,:2].inv()*M[:2,2:]
        sub.update(zip(M[2:,2:],bottom))
    return {k:s.factor(v.xreplace(sub)) for k,v in w.items()}


def quadric(shape,t,b,c,A,B,C,D):
    if shape==76:
        return t**2*B*D+c*C*(A-B)
    return -t**2*B*D+c*A*(B-C)+t*b*(1-c)*B**2+t*c*B*C


def linear_system(shape):
    t,b,c,a = T,BETA,GAMMA,ALPHA
    w = template(shape)
    M = [matrix(w,p,h) for p,h in product(range(2),repeat=2)]
    for m in M:
        assert all(s.factor(block.det()).subs(t,1) in (-1,1) and
                   not any(f.is_Add for f in s.preorder_traversal(s.factor(block.det())))
                   for block in (m[:2,2:],m[2:,:2]))
    bottoms = [scale*m[2:,:2]*m[:2,:2].adjugate()*m[:2,2:] for scale,m in zip((1,a,b,c),M)]
    sub = dict(zip([v for m in M[:2] for v in m[2:,2:]],[v for m in bottoms[:2] for v in m]))
    eq = [f.xreplace(sub)-g for m,B0 in zip(M[2:],bottoms[2:]) for f,g in zip(m[2:,2:],B0)]
    L = s.linear_eq_to_matrix(eq,X[:8])[0]
    factor = 1 if shape==76 else t**-4
    assert s.factor(L.det(method='domain-ge')-factor*(a*t**2-b*c)**4)==0
    L = L.subs(a,b*c/t**2)
    vals = top(shape,t,b,c,A,B,C,D)
    linear = (L*s.Matrix(vals)).applyfunc(s.factor)
    delta = t**2+b*c-b-c
    expected = [0,0,0,0,0,B*delta,0 if shape==76 else B*delta,
                -c*delta*(C if shape==76 else A-2*t*B+t*C)/t**2]
    assert all(s.factor(v-e)==0 for v,e in zip(linear,expected))
    # These four rows solve uniquely for x1,x0,x7,x3 over the parameter torus.
    minor = L.extract((3,1,2,0),(1,0,7,3)).det()
    assert minor!=0 and not any(v.is_Add for v in s.preorder_traversal(s.factor(minor)))
    dets = [s.factor(m[:2,:2].det().xreplace(dict(zip(X,vals)))) for m in M]
    ratios = (b*c/t**2,b,c)
    constraints = [s.factor(dets[0]-r*d) for r,d in zip(ratios,dets[1:])]
    if shape==112:
        off = {B:0,A:-t*C}
        for f,g in zip(constraints,[C**2*(b-c),C**2*(b-1),C**2*(t**2-c)]):
            unit = s.factor(f.subs(off)/g)
            assert unit!=0 and unit.free_symbols<={t,b,c}
            assert not any(v.is_Add for v in s.preorder_traversal(unit))
    else:
        assert dets[0].subs({B:0,C:0})==0
    F = quadric(shape,t,b,c,A,B,C,D)
    for f,h in zip(constraints,(-(b-c),-c*(b-1),b*(c-1))):
        if shape==112:
            h = -h
        assert s.rem(s.fraction(s.factor(f-h*F/t**2))[0],delta,t)==0
    return M


def gauge(shape):
    w = template(shape)
    keys = [k for k in w if w[k] not in X]
    def row(k):
        return [int(k[p]==j) for p,j in product(K,repeat=2)]
    d,_,r = smith(s.Matrix(list(map(row,keys))))
    ds = [n for n in d.diagonal() if n]
    assert ds==([1]*11 if shape==76 else [1]*11+[2])
    holes = {v:k for k,v in w.items() if v in X}
    effects = s.Matrix([row(holes[x]) for x in X[:8]])*r
    if shape==112:
        assert all(v%2==0 for v in effects[:,11])
    groups = ((0,1,2,4),(3,5,6,7)) if shape==76 else (tuple(range(8)),)
    eff = effects[:,len(ds):]
    assert eff.rank()==len(groups)
    assert all(eff.row(i)==eff.row(group[0]) for group in groups for i in group)


def invariants(shape):
    t,b = s.Integer(2),s.Integer(3)
    c = (t**2-b)/(1-b)
    DD = s.solve(quadric(shape,t,b,c,7,1,5,D),D)[0]
    W = completed(shape,t,b,c,7,1,5,DD)
    w = {k:W.get(k,0) for k in WORDS}
    assert all(s.Matrix([[w[i,j,i^h,j^h] for j in K] for i in K]).det()!=0 for h in K)
    def residue(v):
        num,den = s.fraction(v)
        return int(num)*pow(int(den),-1,1009)%1009
    lie = [[residue(w.get(key[:p]+(j,)+key[p+1:],0)) if key[p]==i else 0
            for p,i,j in product(K,repeat=3)] for key in product(K,repeat=4)]
    assert nmod_mat(lie,1009).rank()==61
    G = alternating_gradient(w)
    cap = sum(w[k]*G[k] for k in WORDS)/4
    inv = []
    for p in K:
        L = [matrix(w,p,h)*matrix(G,p,h).T/cap for h in K]
        inv.append(s.factor(sum(sum((m*m)[i,i] for m in L)**2 for i in K)))
    assert len(set(inv))==4


if __name__=='__main__':
    for shape in (76,112):
        linear_system(shape)
        gauge(shape)
        invariants(shape)
    print('PASS: all corner equations, torus stabilizers, generic opposite rank, finite stabilizer and period four')
