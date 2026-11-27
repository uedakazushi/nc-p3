"""Cover the rank-ten surfaces by free presentations and four finite residue fields."""
import json
import os
import subprocess
from itertools import product

from verify_first_rectangle_edge_center import edge_rules
from verify_first_rectangle_edge_presentation import CERTIFICATE
from verify_parallel_curve_residues import residue_test
from verify_parallel_invertible_rectangles import s


def check():
    v,z = s.symbols('v z')
    unit = v*z*(v-1)*(v+1)*(v+z)*(z-1)*(v*z+1)
    line = v+2*z-1
    cert = json.loads(CERTIFICATE.read_text())
    sources = {}
    for key,values in cert.items():
        factors = [s.sympify(value.replace('^','**')) for value in values]
        if int(key[-1])>=2:
            factors = [s.fraction(s.factor(f.subs(v,1/v)))[0] for f in factors]
        sources[key] = factors

    def fmt(f):
        return str(f).replace('**','^')

    for start,keys in enumerate((('A0','B0','A2'),('A3','B1','A1','B3'))):
        a,b,*others = (sources[k] for k in keys)
        # Rabinowitsch eliminates the boundary in every version of Singular.
        code = ['ring r=0,(v,z,T),dp;ideal J;poly q;poly boundary=1-T*('+fmt(unit*(line if start else 1))+');']
        for f,g in product(a,b):
            code.append('J=std(ideal('+fmt(f)+','+fmt(g)+',boundary));')
            for fs in others:
                code.append('if(J[1]!=1){q=1;')
                code.extend('q=reduce(q*('+fmt(h)+'),J);' for h in fs)
                code.append('J=std(J+ideal(q));}')
            code.append('print(J[1]);')
        res = subprocess.run([os.environ.get('SINGULAR_BIN','Singular'),'-q'],input='\n'.join(code+['quit;']),
                             text=True,capture_output=True,check=True)
        assert [r for r in res.stdout.splitlines() if r.strip() and not r.startswith('//')]==['1']*(len(a)*len(b))
    # The reciprocal A3 chart covers the rest of the line except for four number fields.
    on_line = (1-v)/2
    fields = {s.Poly(q,v).monic().as_expr() for f in sources['A3']
              for q,_ in s.factor_list(f.subs(z,on_line))[1] if q!=v+1}
    samples = ((v**3-4*v**2+3*v-4,1013,412),(v**2-2,1009,570),
               (v**2+v+1,1009,634),(2*v**2-3*v+2,1009,274))
    assert fields=={s.Poly(f,v).monic().as_expr() for f,_,_ in samples}
    def modular(value,p):
        a = s.Integer(value)
        def convert(x):
            n,d = s.fraction(x)
            return int(n)*pow(int(d),-1,p)%p
        return edge_rules(a,on_line.subs(v,a),convert,(0,0,0,0))
    for f,p,r in samples:
        assert f.subs(v,r)%p==0 and f.diff(v).subs(v,r)%p!=0
        assert s.gcd(int(s.fraction(unit.subs(z,on_line).subs(v,r))[0]),p)==1
        residue_test(r,p,modular)
    print('PASS: all 35 and 95 chart intersections and all four remaining characteristic-zero factors')


if __name__=='__main__':
    check()
