"""Unmarked generic orbits on the three large horizontal 112 strata."""
from itertools import product

from verify_parallel_invertible_rectangles import K, T, X, bit, candidates, completed, s, template
from verify_third_rectangle_orbits import support_maps


t,b,a,y = s.symbols('t b a y',nonzero=True)
m = t**2
c = (m-b)/(1-b)


def window(t,b,c,a,y):
    d = (c*a*(1-y)+b*(1-c)+c*y)/t
    return completed(112,t,b,c,t*a,1,y,d)


def geometry():
    den = (b-a)*(b-1+a)
    cc = b**2/den
    assert s.factor(b+cc-b*cc-b*a*(1-a)/den)==0
    W2 = b*a*(1-a)*den
    k = -b**2*(1-b)
    x = k/a-k
    A,B = b*(1-b)**2,b**3
    assert s.factor(k**2*W2/a**4-x*(x+A)*(x+B))==0
    c4 = 16*(A**2-A*B+B**2)
    disc = 16*A**2*B**2*(A-B)**2
    assert s.factor(disc-16*b**10*(b-1)**4*(2*b-1)**2)==0
    for v,orders in ((0,(2,10)),(1,(0,4)),(s.Rational(1,2),(0,2))):
        assert tuple(next(n for n in range(11) if s.diff(f,b,n).subs(b,v)!=0)
                     for f in (c4,disc))==orders
    assert (8-s.degree(c4,b),24-s.degree(disc,b))==(2,8)
    norm = m*(a-b)*(y-b)-b*(a-1)*(y-1)
    sol = b*(m*(a-b)-(a-1))/(m*(a-b)-b*(a-1))
    assert s.factor(norm.subs(y,sol))==0
    assert s.gcd(s.expand(norm).coeff(y),norm.subs(y,0))==1


def transformations():
    w = window(t,b,c,a,y)
    holes = [k for k,v in template(112).items() if v in X]
    thin = [k for k in w if k not in holes]
    def weight(k):
        return [int(k[p]==i) for p,i in product(K,repeat=2)]
    def exponent(f):
        e = s.factor(t*s.diff(f,t)/f)
        assert e.is_Integer and f==t**e
        return e
    gauge = s.Matrix([weight(k)+[exponent(w[k])] for k in thin])
    flip = (t,m/b,m/c,1-c+c*y/b,1-c+c*a/b)
    step = (1/t,c/m,1/b,c*(1-a)/m,c*(1-y)/m)
    half = (t,b,c,y,a)
    orig = (t,b,c,a,y)
    def compose(expr,point):
        return tuple(s.factor(f.xreplace(dict(zip((t,b,a,y),(point[0],point[1],point[3],point[4]))))) for f in expr)
    def same(left,right):
        return all(s.factor(z-q)==0 for z,q in zip(left,right))
    assert same(compose(flip,flip),orig)
    assert same(compose(step,step),compose(flip,half))
    assert same(compose(step,flip),compose(flip,step))
    u = (b-t)/(b+t)
    U0,U1 = (a+flip[3])/2,(a-flip[3])/(2*u)
    assert same(compose((u,U0,U1),flip),(-u,U0,U1))
    assert s.factor(b-t*(1+u)/(1-u))==0
    assert s.factor((U0-u*U1-1+c)*b/c-y)==0
    size = {144:(1,0,1,0),146:(1,0,1,0),147:(2,2,2,2)}
    perm_by_shift = []
    for idx,counts in size.items():
        support = [k for k in w if candidates(1)[idx]&bit(k)]
        perms = [list(support_maps([k[j:]+k[:j] for k in support],support)) for j in K]
        assert tuple(map(len,perms))==counts
        if idx!=147:
            assert perms[0][0]==(tuple(K),)*4
            assert perms[2][0]==((1,0,2,3),)*4
        else:
            perm_by_shift = perms

    for shift,destinations in enumerate(((orig,flip),(step,compose(flip,step)),(half,compose(flip,half)))):
        used = set()
        for perm in perm_by_shift[shift]:
            z = {tuple(p[i] for p,i in zip(perm,k[shift:]+k[:shift])):f for k,f in w.items()}
            sol = next(iter(s.linsolve((gauge,s.Matrix([exponent(z[k]) for k in thin])))))
            subs = dict.fromkeys(set().union(*(v.free_symbols for v in sol)),0)
            exps = [v.subs(subs) for v in sol]
            tt = t**exps[-1]
            def raw(k):
                return z[k]*t**(-sum(exps[4*p+i] for p,i in enumerate(k)))
            assert all(s.factor(raw(k)-template(112)[k].subs(T,tt))==0 for k in thin)
            scale = raw(holes[4])
            entries = [s.factor(raw(k)*(1/scale if j<8 else scale)) for j,k in enumerate(holes)]
            coords = (tt,entries[1],(tt**2-entries[1])/(1-entries[1]),entries[2]/tt,entries[5])
            j, = [j for j,point in enumerate(destinations) if same(coords,point)]
            used.add(j)
            target = window(*destinations[j])
            assert all(s.factor(expr-target[k])==0 for k,expr in zip(holes,entries))
        assert used=={0,1}


if __name__=='__main__':
    geometry()
    transformations()
    print('PASS: elliptic K3 invariants, rational slice, invariant field and all window transformations')
