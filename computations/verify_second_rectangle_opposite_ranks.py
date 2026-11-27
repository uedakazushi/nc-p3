"""Opposite ranks on all large horizontal 112 strata and the rank-fourteen elliptic surface."""
from itertools import combinations, product

from verify_parallel_invertible_rectangles import K, s
from verify_parallel_horizontal_moduli import window


def opposite_minors():
    u,v,z,w = s.symbols('u v z w')
    t = (u+v)/(1+u*v)
    b,c = u*t,v*t
    x,y = u*(z+v)/(1+u*v),u*(w+v)/(1+u*v)
    tensor = window(t,b,c,x,y)
    unit = u*v*(u+v)*(1+u*v)*(z*w-1)*(z+v)*(w+v)*(z-u)*(w-u)*(u*z-1)*(u*w-1)*(v*z+1)*(v*w+1)
    variables = (u,v,z,w)
    def monic(f):
        return s.Poly(f,*variables).monic().as_expr()
    units = {monic(f) for f,_ in s.factor_list(unit)[1]}
    def strip(value):
        assert value!=0
        a,d = s.fraction(s.factor(value))
        assert all(monic(f) in units for f,_ in s.factor_list(d)[1])
        return s.prod((f**e for f,e in s.factor_list(a)[1] if monic(f) not in units),start=s.S.One)
    C = (u+v)*(1+u*v)
    H,G = u*(v+1)**2,v*(u+1)**2
    F = C*(z+w-2)-H*(z*w-1)
    L = C*(z+w+2)+G*(z*w-1)
    for h,(divisor,extra) in enumerate(((v-1,F),(u-1,L),((u-1)*(v-1),1),((u-1)*(v-1),1))):
        M = s.Matrix([[tensor.get((i,j,i^h,j^h),0) for j in K] for i in K])
        assert s.factor(strip(M.det(method='domain-ge'))/(divisor**2*extra)) in (-1,1)
        for n,base in ((2,1),(3,divisor)):
            polys = [s.factor(strip(d)/base) for I,J in product(combinations(K,n),repeat=2)
                     if (d:=s.factor(M.extract(I,J).det(method='domain-ge')))!=0]
            assert all(s.fraction(p)[1]==1 for p in polys)
            if n==2 or h>=2:
                assert 1 in polys or -1 in polys,(h,n,polys)
            else:
                assert all(a in polys or -a in polys for a in ((w-1,z-1) if h==0 else (w+1,z+1)))
    assert s.factor(F.subs(u,1)+(1+v)**2*(z-1)*(w-1))==0
    assert s.factor(L.subs(v,1)-(1+u)**2*(z+1)*(w+1))==0
    for f in (F,L):
        f = s.Poly(f,w)
        assert s.gcd(*f.all_coeffs())==1
    N = G+H
    P,Q = 2*(G-H)/N,1-4*C/N
    A,B = u+1/u,v+1/v
    delta = A**2+A*B+B**2+2*A+2*B-4
    assert s.factor(P**2-4*Q-16*(u*v)**2*delta/N**2)==0


def elliptic_surface():
    A,v = s.symbols('A v')
    b = A+2
    c = A**2+2*A-2
    model = v**4+b*v**3+c*v**2+b*v+1
    I = 12-3*b**2+c**2
    J = 72*c+9*b**2*c-54*b**2-2*c**3
    disc = (A-2)**3*(A+2)**3*(3*A+10)**2
    assert s.factor(s.discriminant(model,v)/disc)==1 and s.factor((4*I**3-J**2)/disc)==27
    assert s.degree(I,A)==4 and s.degree(J,A)==6
    assert s.rem(I,A-2,A)==0 and s.rem(I,(A-2)**2,A)!=0
    assert all(I.subs(A,x)!=0 for x in (-2,-s.Rational(10,3)))
    u = s.Symbol('u')
    assert s.factor(3*(u+1/u)+10-(3*u+1)*(u+3)/u)==0
    assert s.factor(((u+1/u)**2-4)-(u-1/u)**2)==0
    assert 2*3*2*4==48 and 4*6*2**2*4**2==1536


def period_two_curve():
    l,m = s.symbols('l m')
    a,b = 2*(1+l*l)/(1-l*l),2*(1+m*m)/(1-m*m)
    d = s.factor((a*a+a*b+b*b+2*a+2*b-4)*(1-l*l)**2*(1-m*m)**2/16)
    assert s.expand(d-1+3*l*l*m*m-l*l*m*m*(l*l+m*m))==0
    r = l*m
    S,E = r*(l+m)/(r+1),r*(l-m)/(r-1)
    for expr in (S*S-(2*r-1),E*E+2*r+1,(a-b)/(a+b+4)+S*E):
        assert s.rem(s.fraction(s.factor(expr))[0],d,l,m)==0
    assert s.factor(((r+1)*S+(r-1)*E)/(2*r)-l)==0
    assert s.factor(((r+1)*S-(r-1)*E)/(2*r)-m)==0
    # Exceptional inverse denominators supply no additional component.
    for f in (r,r-1,r+1,l-1,l+1,m-1,m+1):
        assert s.gcd(f,d)==1


if __name__=='__main__':
    opposite_minors()
    elliptic_surface()
    period_two_curve()
    print('PASS: all opposite-rank minors, elliptic surface discriminant, and rational period-two curve')
