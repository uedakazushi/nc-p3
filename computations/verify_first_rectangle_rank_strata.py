"""All opposite ranks and the intersections of the support-85 rank loci."""
from itertools import combinations,product

from verify_parallel_invertible_rectangles import K,completed,s


def check():
    u,v,z,r = s.symbols('u v z r')
    C = (u+v)*(1+u*v)
    t = (u+v)/(1+u*v)
    b,c = u*t,v*t
    x = u*(z+v)/(1+u*v)
    w = completed(76,t,b,c,x,1,1,v*(1-u*z)/(u+v))
    units = u*v*z*C*(u*z-1)*(z-u)*(v*z+1)*(v+z)
    F,G = u*(v+1)**2*z-C,v*(u+1)**2*z+C
    def strip(poly):
        f,_ = s.fraction(s.factor(poly))
        return s.prod(g**e for g,e in s.factor_list(f)[1] if s.div(units,g,u,v,z)[1]!=0)
    tests = [(v-1,F),(u-1,G),((u-1)*(v-1),1),((u-1)*(v-1),1)]
    for h,(L,Q) in enumerate(tests):
        M = s.Matrix([[w.get((i,j,i^h,j^h),0) for j in K] for i in K])
        det = strip(M.det())
        assert s.factor(det/(L**2*Q)) in (-1,1)
        for n,g in ((3,L),(2,1)):
            fs = [strip(f) for I,J in product(combinations(K,n),repeat=2) if (f:=M.extract(I,J).det())!=0]
            if g!=1:
                assert all(s.rem(f,g,u,v,z)==0 for f in fs)
            assert any(s.factor(f-g)==0 or s.factor(f+g)==0 for f in fs)
    assert s.factor(F.subs(u,1)-(v+1)**2*(z-1))==0
    assert s.factor(G.subs(v,1)-(u+1)**2*(z+1))==0
    U,V = r*(r-1)/(r+1),(1-r)/(r*(r+1))
    Z = -4*V/(V+1)**2
    sub = {u:U,v:V,z:Z}
    assert all(s.factor(f.xreplace(sub))==0 for f in (F,G))
    assert s.factor(units.xreplace(sub))!=0 and s.factor((U+1)/(V+1)-r)==0
    # The positive t=1 strata of both horizontal shapes are one-sided symmetric.
    a,y = s.symbols('a y')
    for b,c in ((s.S.One,v),(u,s.S.One)):
        for shape in (76,112):
            W = (completed(76,1,b,c,a,1,1,c*(1-a)) if shape==76 else
                 completed(112,1,b,c,a,1,y,c*a*(1-y)+b*(1-c)+c*y))
            for (i,j,k,l),f in W.items():
                target = (k,j,i,l) if b==1 else (i,l^1,k,j^1)
                assert s.factor(f-W.get(target,0))==0
    print('PASS: all opposite rank minors, rational rank-14 intersection and both one-sided symmetries')


if __name__=='__main__':
    check()
