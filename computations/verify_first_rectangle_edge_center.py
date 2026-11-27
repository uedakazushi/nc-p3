"""A finite module over the five quartic central forms at every rank-ten 76 edge point."""
from itertools import product
import os
import subprocess

from flint_rational_function import rational_function_field
from verify_parallel_invertible_rectangles import K,completed,matrix,s
from verify_rectangular_surface_cut_charts import chart,relation_cut_data


def edge_rules(v,z,convert,shifts):
    x = (z+v)/(1+v)
    w = completed(76,1,1,v,x,1,1,v*(1-x))
    rel = []
    for p in K:
        q = (p+1)%4
        rows = []
        for h in K:
            for col in matrix(w,p,h).columnspace():
                rows.append((h^shifts[p]^shifts[q],{(i^shifts[p],i^h^shifts[q]):convert(s.factor(col[i]))
                                                     for i in K if col[i]}))
        rel.append(rows)
    return rel


def check():
    v,z = s.symbols('v z')
    units = v*z*(v-1)*(v+1)*(v+z)*(z-1)*(v*z+1)
    F = rational_function_field(('v','z'))
    rel = edge_rules(v,z,F.from_sympy,(0,2,2,0))
    data = relation_cut_data(rel,F,F.from_sympy(units).n,full_center=True)
    profile = (1,4,10,20,30,36,34,20,1,0)
    factors = [chart(p,data,profile) for p in K]
    assert list(map(len,factors))==[3,16,24,3]
    A,B,C,D = [[s.sympify(f.replace('^','**')) for f in fs] for fs in factors]
    def reciprocal(fs):
        return [s.fraction(s.factor(f.subs(v,1/v)))[0] for f in fs]
    def fmt(f):
        return str(f).replace('**','^')
    code = ['LIB "elim.lib";ring r=0,(v,z),dp;ideal J;poly u='+fmt(units)+';']
    for a,d in product(A,D):
        code.append(f'J=std(sat(ideal({fmt(a)},{fmt(d)}),ideal(u))[1]);')
        for fs in (reciprocal(A),reciprocal(D),B,C,reciprocal(B),reciprocal(C)):
            g = '*'.join('('+fmt(f)+')' for f in fs)
            code.append(f'if(J[1]!=1){{J=std(J+ideal(reduce({g},J)));}}')
        code.append('J=std(sat(J,ideal(u))[1]);print(J);')
    res = subprocess.run([os.environ.get('SINGULAR_BIN','Singular'),'-q'],input='\n'.join(code+['quit;']),
                         text=True,capture_output=True,check=True)
    assert [v for v in res.stdout.splitlines() if v.strip() and not v.startswith('//')]==['1']*9
    print('PASS: five central forms annihilate all length-twelve paths in the rank-ten edge quotient')


if __name__=='__main__':
    check()
