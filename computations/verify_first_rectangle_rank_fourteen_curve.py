"""Certify every point of the opposite-rank-fourteen curve on horizontal support 85."""
from itertools import combinations
import json
from pathlib import Path

from flint_rational_function import rational_function_field
from localized_unit_echelon import unit_echelon
from quartic_center_presentation import certify_presentation
from verify_parallel_curve_residues import residue_test
from verify_parallel_invertible_rectangles import K,completed,matrix,s
from verify_rectangular_surface_cut_charts import relation_cut_data


def curve(r):
    u = r*(r-1)/(r+1)
    v = (1-r)/(r*(r+1))
    z = -4*v/(v+1)**2
    t = s.factor((u+v)/(1+u*v))
    b,c = s.factor(u*t),s.factor(v*t)
    x,d = s.factor(u*(z+v)/(1+u*v)),s.factor(v*(1-u*z)/(u+v))
    return completed(76,t,b,c,x,1,1,d),(u+v,1+u*v,u-1,v-1)


def rules(w,shifts,convert):
    return [[(h^shifts[j]^shifts[(j+1)%4],
              {(i^shifts[j],i^h^shifts[(j+1)%4]):convert(s.factor(col[i])) for i in K if col[i]})
             for h in K for col in matrix(w,j,h).columnspace()] for j in K]


def check():
    r = s.Symbol('r')
    w,extra = curve(r)
    field = rational_function_field(('r',))
    unit = s.prod({f for e in tuple(w.values())+extra for q in s.fraction(s.factor(e)) for f,_ in s.factor_list(q)[1]})
    units = field.from_sympy(unit).n
    for j in K:
        for h in K:
            M = s.Matrix.hstack(*matrix(w,j,h).columnspace())
            numerators = [s.fraction(s.factor(M.extract(J,range(M.cols)).det()))[0] for J in combinations(K,M.cols)]
            gcd = s.gcd_list(numerators)
            assert all(s.rem(unit,f,r)==0 for f,_ in s.factor_list(gcd)[1])
    charts = []
    for cc,shifts in enumerate(((0,2,2,0),(2,0,2,0))):
        factors = set()
        def record(f):
            factors.update(str(p) for p,_ in f.factor()[1])
        def span(rows):
            return unit_echelon(rows,units,record)
        def pivot(rows):
            indices = []
            return unit_echelon(rows,units,record,indices),indices
        rel = rules(w,shifts,field.from_sympy)
        _,_,_,sections = relation_cut_data(rel,field,units,True,span)
        certify_presentation(rel,sections,0,[field.one*int(i==cc) for i in range(5)],field.one,
                             pivot,span,(1,4,10,20,30,36,34,20,0))
        charts.append(factors)
    assert list(map(len,charts))==[120,141]
    expected = json.loads((Path(__file__).resolve().parents[1]/'research'/'first_rectangle_rank_fourteen_residues.json').read_text())
    assert set.intersection(*charts)==set(expected)
    for f,(p,a) in expected.items():
        poly = s.Poly(s.sympify(f.replace('^','**')),r)
        assert poly.is_irreducible and poly.eval(a)%p==0 and poly.diff().eval(a)%p and unit.subs(r,a)%p
        def modular(value,prime):
            W,_ = curve(s.Integer(value))
            def convert(x):
                n,d = s.fraction(x)
                return int(n)*pow(int(d),-1,prime)%prime
            return rules(W,(0,0,0,0),convert)
        residue_test(a,p,modular)
    print('PASS: both presentations and all eight number fields of the opposite-rank-fourteen curve')


if __name__=='__main__':
    check()
