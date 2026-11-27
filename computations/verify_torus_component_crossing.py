"""The three fixed directions and the full tangent at the meeting of two torus families."""
import sympy as sp

from verify_finite_torus_slices import normal_table, evaluation_rank, simple_local
from verify_simple_weight_torus import A,B,C,D,H,F,X,Y,monomial_chart,schur,solve_chart


def verify():
    chart,_,_ = monomial_chart()
    main,_ = solve_chart(chart)
    w = {q:v.subs({B:sp.Rational(1,2),C:1}) for q,v in main.items()}
    assert normal_table(w,(0,1,2,3),6)==[0,1,2,2,1,1,3,1,1,2,2,1,0]
    assert evaluation_rank(lambda r,i:simple_local(w,r,i))==19
    eq = schur(chart)
    eq = [p.subs({X:H/C+1/A,Y:B*D/C+D*H/(A*B)}) for p in eq]
    h = B*C**2*(1-B)/(C-B)
    sub = {D:A*(C-B)/C,H:h,F:1+A*(h-B)/C}
    u,v,z = sp.symbols('u v z')
    c = 1+z
    a = u+v+2*c
    b = (u+c)/a
    reduced = [sp.factor(p.subs(sub,simultaneous=True).subs({A:a,B:b,C:c},simultaneous=True)) for p in eq]
    zero = {u:0,v:0,z:0}
    assert all(sp.fraction(p)[1].subs(zero)!=0 for p in reduced)
    assert sp.factor(reduced[1]/(u*v)).subs(zero)!=0
    for p in reduced:
        assert sp.cancel(p.subs(u,0))==0 and sp.cancel(p.subs({v:0,z:0}))==0
    quotient = sp.factor(reduced[6].subs(v,0)/(u*z))
    num,den = sp.fraction(quotient)
    assert num.subs(zero)!=0 and den.subs(zero)!=0


if __name__ == '__main__':
    verify()
    print('The fixed slice ideal is (uv,uz), with 17 unrestricted transverse tangent dimensions and 19 identities.')
