"""A common torus chart after removal of the anti-diagonal and all low 76 fibers."""
from itertools import combinations,product

from flint_rational_function import rational_function_field
from localized_unit_echelon import unit_echelon
from verify_klein_center_lift import central_equations,field_paths
from verify_parallel_invertible_rectangles import K,completed,matrix,s


def coordinates():
    u,v,z,w = s.symbols('u v z w')
    t = (u+v)/(1+u*v)
    b,c = u*t,v*t
    x,y = u*(z+v)/(1+u*v),u*(w+v)/(1+u*v)
    E = (b+c)*x-b*c
    assert s.factor(t*t-b-c+b*c)==0
    assert s.factor(E-b*t*z)==0 and s.factor(((b+c)*y-b*c)-b*t*w)==0
    assert s.factor((b+c)*x*y-b*c*(x+y)+b*b*(c-1)-u*u*(u+v)**2*(z*w-1)/(1+u*v)**3)==0
    def evaluate(expr,point):
        return tuple(s.factor(f.xreplace(dict(zip((u,v,z,w),point)))) for f in expr)
    param = (t,b,c,x,y)
    for action,dest in (((v,1/u,-z,-w),(1/t,c/t**2,1/b,c*(1-x)/t**2,c*(1-y)/t**2)),
                        ((1/u,1/v,w,z),(t,t**2/b,t**2/c,1-c+c*y/b,1-c+c*x/b))):
        assert all(s.factor(f-g)==0 for f,g in zip(evaluate(param,action),dest))
    raw = completed(76,t,b,c,x,1,1,v*(1-u*z)/(u+v))
    factors = {f for a in raw.values() for p in s.fraction(s.factor(a)) for f,_ in s.factor_list(p)[1]}
    expected = (u,v,z,u+v,1+u*v,z+v,1-u*z,z-u,1+v*z)
    assert {s.Poly(f,u,v,z).monic().as_expr() for f in factors}=={s.Poly(f,u,v,z).monic().as_expr() for f in expected}


def low_window(t,b,c,x):
    raw = completed(76,t,b,c,x,1,1,c*(1-x)/t**2)
    variables = tuple(sorted(set().union(*(v.free_symbols for v in raw.values())),key=str))
    field = rational_function_field(tuple(map(str,variables)))
    units = s.prod({f for z in list(raw.values())+[t,b,c,b-1,b+c] for q in s.fraction(s.factor(z)) for f,_ in s.factor_list(q)[1]})
    U = field.from_sympy(units).n
    def unit(f):
        for poly in s.fraction(s.factor(f)):
            if not poly or any(not s.div(units,g,*variables)[1]==0 for g,_ in s.factor_list(poly)[1]):
                return False
        return True
    rel = []
    for p in K:
        rules = []
        for h in K:
            columns = matrix(raw,p,h).columnspace()
            frame = s.Matrix.hstack(*columns)
            assert any(unit(frame.extract(I,range(len(columns))).det()) for I in combinations(K,len(columns)))
            for col in columns:
                rules.append((h,{(j,j^h):field.from_sympy(s.factor(col[j])) for j in K if col[j]}))
        rel.append(rules)
    def echelon(rows):
        return unit_echelon(rows,U)
    basis,image = field_paths(rel,field.one,echelon)
    profile = (1,4,10,20,35,56)
    assert all(sum(len(basis[p,n,h]) for h in K)==profile[n] for p,n in product(K,range(6)))
    for h,d in enumerate((3,2,0,0)):
        rows,n = central_equations(basis,image,h)
        assert n-len(echelon(rows))==d


if __name__=='__main__':
    coordinates()
    t,b,c,x = s.symbols('t b c x')
    low_window(t,b,(t*t-b)/(1-b),x)
    for sign in (s.S.One,-s.S.One):
        low_window(sign,s.S.One,c,x)
    print('PASS: common torus coordinates, rotations and low paths and quartic centers at every candidate on 85')
