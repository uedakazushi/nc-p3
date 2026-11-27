"""The two intrinsic lines, an indexed parameter and the shift involution on the PBW curve."""
from itertools import product

import sympy as sp
from verify_klein_mixed_supports import smith, character
from verify_klein_single_missing_windows import regular_window, relations, K


def fact(m):
    return m.applyfunc(sp.factor)


def transport(w,p):
    def split(q):
        return ((0,1),(2,3)) if q%4<2 else ((2,3),(0,1))
    a,b = split(p)
    c,d = split(p+1)
    rel = relations(w,p)
    def block(r,left,right):
        return sp.Matrix([[r.get((i,j),0) for j in right] for i in left])
    j = sp.Matrix([[0,-1],[1,0]])
    f = j*next(block(r,a,c) for r in rel if set(r)&set(product(a,c))).T
    g = j*next(block(r,b,d) for r in rel if set(r)&set(product(b,d))).T
    mixed = [r for r in rel if set(r)&set(product(b,c))]
    left = sp.Matrix([list((block(r,b,c)*f).T) for r in mixed])
    right = sp.Matrix([list(block(r,a,d)*g) for r in mixed])
    return f,g,fact(left.inv()*right)


def unsigned_monomials(w,t):
    first,second,third,fourth = (2,3,0,1),(3,2,1,0),(0,1,3,2),(1,0,2,3)
    common = (2,3,1,0)
    opposite = (3,2,0,1)
    keys = sorted(w)
    exponents = sp.Matrix([[int(k[p]==i) for p,i in product(K,repeat=2)] for k in keys])
    d,l,_ = smith(exponents)
    rk = sum(x!=0 for x in d.diagonal())
    for parameter,rotation,reverse,perms in [
        (-t,0,False,(first,first,second,second)),(1/t,0,False,(third,third,fourth,fourth)),
        (t,2,False,(common,)*4),(t,0,True,(common,common,opposite,opposite)),
    ]:
        target = {}
        for k,v in w.items():
            word = k[::-1] if reverse else k
            target[word[rotation:]+word[:rotation]] = v.subs(t,parameter)
        ratio = [w[k]/target[tuple(pi[i] for pi,i in zip(perms,k))] for k in keys]
        assert all(character(row,ratio)==1 for row in l.tolist()[rk:])


def main():
    t,z = sp.symbols('t z',nonzero=True)
    w = regular_window(t)
    unsigned_monomials(w,t)
    L = 1+t**2
    H = t+1/t
    X = H**2
    form = sp.Matrix([[0,0,0,1],[0,0,-1,0],[0,-1,0,0],[1,0,0,0]])
    f,g = sp.eye(2),sp.eye(2)
    normal = []
    for p in K:
        a,b,m = transport(w,p)
        k = fact(m*m)[0,0]
        assert fact(m*m)==k*sp.eye(4) and sp.factor(m.trace())==0
        for plus in (True,False):
            plane = sp.Matrix([[1,0],[0,1],[z,0],[0,z]]) if plus else sp.Matrix([[1,0],[z,0],[0,1],[0,z]])
            quartic = sp.Poly(sp.factor(plane.row_join(m*plane).det()),z)
            assert quartic.degree()==4 and quartic.nth(3)==quartic.nth(1)==0
            normalized = quartic.nth(2)**2/(4*quartic.nth(0)*quartic.nth(4))
            assert sp.factor(normalized-(H+(1 if plus else -1))**2)==0
            for v in (quartic.nth(0),quartic.nth(4),k):
                for poly in sp.fraction(sp.factor(v)):
                    assert all(c in (t,L) for c,_ in sp.factor_list(poly)[1])
        adj = form*m.T*form
        invariant = sp.factor(((adj*m/k)**2).trace())
        assert sp.factor(invariant-(4-32/X if p%2==0 else 4-2*X))==0
        change = sp.kronecker_product(f,g)
        normal.append(fact(change.inv()*m*change))
        f,g = fact(a*f),fact(b*g)
    scale = -4*t**3/L
    assert fact(normal[2]-scale*normal[0])==sp.zeros(4)
    s = sp.I*(t-sp.I)/(t+sp.I)
    a = sp.Matrix([[1,-sp.I*t/L],[-4*sp.I*t/L,4*t**2/L**2]])
    b = sp.Matrix([[1,-1/t],[-s,-s/t]])
    c = sp.kronecker_product(a,b)
    q = 8*(1-sp.I)*t**2*(t-sp.I)/(t+sp.I)**4
    r = (1+sp.I)/(2*t*(t-sp.I)**2*(t+sp.I))
    for left,right,k in ((normal[0],normal[1],q),(normal[1],normal[0],r)):
        diff = k*c*left-right.subs(t,s)*c
        assert diff.applyfunc(lambda v:sp.cancel(v,extension=sp.I))==sp.zeros(4)
    assert sp.cancel((s+1/s)**2-16/X,extension=sp.I)==0
    print('PASS: eight endpoint quartics, four phase invariants, fixed-phase identifications, reversal and shift conjugacy')


if __name__=='__main__':
    main()
