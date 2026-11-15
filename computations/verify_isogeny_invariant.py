#!/usr/bin/env python3
"""Exact algebraic checks for the intrinsic elliptic 2-isogeny of nc-p3.

Checks homogeneous cubic quotient identity, both Hesse j functions, the
rational Hauptmodul coordinate and ramification data. Independent of any
numerical sampling. Mathematical intrinsicness is proved in the manuscript.
"""
import argparse
import json
from pathlib import Path
import sympy as sp


def verify():
    w,y,z,t,h,x = sp.symbols('w y z t h x')
    source = w**3+h*y**3+t*z**3+(2-t*h)*w*y*z
    p=sp.Matrix([w,y,z])
    K=sp.Matrix([[y,w,-t*z],[z,-h*y,w],[-w,z,y]])
    q=sp.Matrix([w**2-t*h*y*z, -w*y-t*z**2, -w*z-h*y**2])
    # Cross product of the first two rows of K, so K(p)q=0 on E.
    assert sp.simplify(q-K.row(0).T.cross(K.row(1).T)) == sp.zeros(3,1)
    assert all(sp.expand(v)==0 for v in K*q+sp.Matrix([0,0,source]))
    omega=p.cross(q)
    A=t*z**3-h*y**3
    B=2*z*w**2+h*y**2*w-t*h*y*z**2
    C=-2*y*w**2-t*z**2*w+t*h*y**2*z
    assert all(sp.expand(omega[i]-[A,B,C][i])==0 for i in range(3))
    quotient=t*h*A**3+t*B**3+h*C**3-(t*h+4)*A*B*C
    factor=A*(h*h*t*y**3+h*t*t*z**3-6*h*t*w*y*z+8*w**3)
    assert sp.expand(quotient-source*factor)==0

    # Source E_x and quotient E'_x, in Hesse normal coordinates.
    k_source=-(2-x)**3/(27*x)
    k_quotient=(x+4)**3/(27*x**2)
    def J(k): return sp.factor(27*k*(k+8)**3/(k-1)**3)
    j=J(k_source); jq=J(k_quotient)
    j_expected=(x-2)**3*(x**3-6*x**2+228*x-8)**3/(x*(x-8)**3*(x+1)**6)
    jq_expected=(x+4)**3*(x**3+228*x**2+48*x+64)**3/(x**2*(x-8)**6*(x+1)**3)
    assert sp.cancel(j-j_expected)==0
    assert sp.cancel(jq-jq_expected)==0
    U=4096*(x+1)**3/(x*(x-8)**3)
    assert sp.cancel(j-(U+256)**3/U**2)==0
    assert sp.cancel(jq-(U+16)**3/U)==0
    n,d=sp.fraction(sp.cancel(U))
    assert max(sp.degree(n,x),sp.degree(d,x))==4
    assert sp.gcd(n,d)==1
    dU=sp.factor(sp.diff(U,x))
    assert sp.cancel(dU+4096*(x+1)**2*(x*x+20*x-8)/(x*x*(x-8)**4))==0
    roots=sp.solve(x*x+20*x-8,x)
    assert len(roots)==2 and all(sp.simplify(U.subs(x,r)+64)==0 for r in roots)
    assert sp.simplify(U.subs(x,-1))==0
    assert sp.factor(U+64) == 64*(x*x+20*x-8)**2/(x*(x-8)**3)
    return {
      'result':'PASS',
      'method':'exact QQ[t,h,w,y,z,x] polynomial and rational operations',
      'source_cubic':str(source),
      'quotient_cubic':'t*h*A**3 + t*B**3 + h*C**3 - (t*h+4)*A*B*C',
      'polynomial_quotient_identity': 'G(A,B,C) = F(w,y,z) * A * (h**2*t*y**3 + h*t**2*z**3 - 6*h*t*w*y*z + 8*w**3)',
      'j_source':str(sp.factor(j)),
      'j_quotient':str(sp.factor(jq)),
      'isogeny_modular_parameter':str(sp.factor(U)),
      'modular_formulas': ['j=(U+256)**3/U**2','j_quotient=(U+16)**3/U'],
      'map_degree':4,
      'branch_points':{'U=0':'x=-1 (ramification 3), x=infinity (unramified)',
                       'U=infinity':'x=8 (ramification 3), x=0 (unramified)',
                       'U=-64':'x=-10+6sqrt(3), -10-6sqrt(3) (each ramification 2)'},
      'ramification_sum':6,
      'geometric_meaning':'A4 monodromy for the degree-four x->U cover (manuscript proof)',
      'limitations':['Pair of elliptic invariants is only a generic necessary condition for indexed isomorphism.',
          'Not a classification of the full nine-parameter family or its moduli component.']
    }

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    obj=verify()
    output=json.dumps(obj,ensure_ascii=False,indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(output,encoding='utf-8')
    print(output)
