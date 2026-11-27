"""A sharp fourth-degree obstruction shared by all five large horizontal 112 supports."""
from itertools import combinations,product
import os
import subprocess

from flint_rational_function import rational_function_field
from localized_unit_echelon import unit_echelon
from verify_klein_center_lift import field_paths
from verify_parallel_invertible_rectangles import K,X,completed,matrix,s,template
from verify_third_rectangle_frames import SINGULAR
from verify_third_rectangle_orbits import support_maps

b,t,x,y = s.symbols('b t x y')
EXTRA = b*b*t*t-b*b*x*y+b*b*x+b*b*y*y-b*b+b*t*t*x*y-b*t*t*x-2*b*t*t*y-b*y*y+b*y+t*t*y


def paths():
    c = (t*t-b)/(1-b)
    E = (b+c)*x*y-b*c*(x+y)+b*b*(c-1)
    w = completed(112,t,b,c,t*x,1,y,(c*x*(1-y)+b*(1-c)+c*y)/t)
    variables = (b,t,x,y)
    field = rational_function_field(tuple(map(str,variables)))
    terms = (t,b,c,b-1,b+c,E,x,y,x-b,y-b,x-1,y-1,c*x+b*(1-c),c*y+b*(1-c))
    units = s.prod({f for value in terms for g in s.fraction(s.factor(value)) for f,_ in s.factor_list(g)[1]})
    U = field.from_sympy(units).n
    rel = []
    for p in K:
        rows = []
        for h in K:
            cols = matrix(w,p,h).columnspace()
            frame = s.Matrix.hstack(*cols)
            minors = [s.factor(frame.extract(I,range(len(cols))).det()) for I in combinations(K,len(cols))]
            assert any(m and all(s.div(units,f,*variables)[1]==0 for g in s.fraction(m)
                                 for f,_ in s.factor_list(g)[1]) for m in minors)
            rows.extend((h,{(j,j^h):field.from_sympy(s.factor(col[j])) for j in K if col[j]}) for col in cols)
        rel.append(rows)
    extra = set()
    def record(f):
        extra.update(str(p) for p,_ in f.factor()[1])
    basis,_ = field_paths(rel,field.one,lambda rows:unit_echelon(rows,U,record))
    assert extra=={str(field.from_sympy(EXTRA).n)}
    assert all(sum(len(basis[p,n,h]) for h in K)==(1,4,10,20,35,56)[n] for p,n in product(K,range(6)))


def cover():
    base = template(112)
    B = next(k for k,f in base.items() if f==X[4])
    for r in K:
        for pi in support_maps([k[r:]+k[:r] for k in base],base):
            moved = {tuple(p[i] for p,i in zip(pi,k[r:]+k[:r])):f for k,f in base.items()}
            assert moved[B] in (X[1],X[4],X[13],X[15])
    u,v,z,w = s.symbols('u v z w')
    tt = (u+v)/(1+u*v)
    torus = {t:tt,b:u*tt,x:u*(z+v)/(1+u*v),y:u*(w+v)/(1+u*v)}
    factor = s.factor(EXTRA.xreplace(torus))
    common = u*v*(u+v)*(1+u*v)*(z*w-1)*(z+v)*(w+v)*(z-u)*(w-u)*(u*z-1)*(u*w-1)*(v*z+1)*(v*w+1)
    orbit = []
    first = (u,v,z,w)
    for _ in K:
        U,V,Z,W = first
        orbit.extend((first,(1/U,1/V,W,Z)))
        first = (V,1/U,-Z,-W)
    for sub,variables,extra in (({},(u,v,z,w),(u*u-1)*(v*v-1)),({u:1},(v,z,w),v*v-1),({u:-1},(v,z,w),v*v-1)):
        unit = s.prod(f for f,_ in s.factor_list(s.factor(common.xreplace(sub)*extra))[1])
        tests = []
        for point in orbit:
            if s.factor((point[0]**2-1).xreplace(sub))==0:
                continue
            n,d = s.fraction(s.factor(factor.xreplace(dict(zip((u,v,z,w),point))).xreplace(sub)))
            assert all(s.div(unit,g,*variables)[1]==0 for g,_ in s.factor_list(d)[1])
            tests.append(s.prod(f for f,_ in s.factor_list(n)[1] if s.div(unit,f,*variables)[1]!=0))
        def fmt(f):
            return str(f).replace('**','^')
        code = SINGULAR.replace('poly u=','poly allowed=').replace('gcd(f,u)','gcd(f,allowed)')
        code = code.replace('v,h,b',','.join(map(str,variables))).replace('v*h*(h-1)*b*(b-1)',fmt(unit))
        code = ('LIB "elim.lib";\n'+code+'ideal J=pairIdeal(ideal('+','.join(map(fmt,tests))+'));'
                +'J=std(sat(J,ideal(allowed))[1]);print(J);quit;')
        res = subprocess.run([os.environ.get('SINGULAR_BIN','Singular'),'-q'],input=code,text=True,capture_output=True,check=True)
        assert [z for z in res.stdout.splitlines() if z.strip() and not z.startswith('//')]==['1'],res.stdout


if __name__=='__main__':
    paths()
    cover()
    print('PASS: exactly one path-chart divisor and all orbit saturations on horizontal 112 supports 143 through 147')
