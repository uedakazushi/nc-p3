#!/usr/bin/env python3
"""Exact point-scheme support for the 2-periodic non-Ore family.

The group-theoretic step excluding Möbius symmetries of the generic quintic
uses irreducibility mod 2, factor type (2,3) mod 5, and PGL_2 stabilizers.
See the manuscript for the proof; this program certifies the algebraic input.
"""
from itertools import combinations
import sympy as sp
import json

w,y,z,x,a,k,v,b,s=sp.symbols('w y z x a k v b s')

def matrix(e):
    return sp.Matrix([[y,w,0,-e*a*z],[z,0,w,0],[-k*w,z,y,0],
       [x,0,-v*z,w],[-b*w,x,-s*z,y],[0,0,x,z]])

def verify():
    results={}
    for e in (1,-1):
        M=matrix(e)
        minors={tuple(c):sp.factor(M.extract(c,range(4)).det())
            for c in combinations(range(6),4)}
        assert len(minors)==15
        # Generic isolated points w=1, z != 0. All rank minors vanish on Q_e(z)=0.
        q=e*a*v**2*z**5+4*s*z**2-4*b
        X=-v*z**2/2
        Y=-k/(2*z)-e*a*v*z**3/4
        for index,polynomial in minors.items():
            den=sp.together(polynomial.subs({w:1,x:X,y:Y}))
            numerator=sp.fraction(den)[0]
            assert sp.rem(numerator,q,z)==0,(e,index,sp.factor(sp.rem(numerator,q,z)))
            # The whole line w=z=0 belongs to the point scheme.
            assert sp.expand(polynomial.subs({w:0,z:0}))==0
        # Necessity from selected minors, recorded for the proof.
        assert sp.factor(minors[(1,2,3,5)])==-z**2*(v*z**2+2*w*x)
        assert sp.expand(minors[(0,1,2,5)]-z*(e*a*x*z**2-k*w**3-2*w*y*z))==0
        assert sp.expand(minors[(0,1,4,5)]-z*(e*a*x**2*z-b*w**3+s*w*z**2))==0
        assert sp.factor(minors[(0,1,2,3)].subs(z,0))==-k*w**4
        results[str(e)]={'rank_four_minors':len(minors), 'all_isolated_point_minors_vanish_mod_quintic':True,
              'line_w=z=0_in_point_scheme':True,'projected_quintic':str(q)}
    f=z**5+z**2+1
    mod2=sp.factor_list(sp.Poly(f,z,modulus=2))[1]
    mod5=sp.factor_list(sp.Poly(f,z,modulus=5))[1]
    assert [(p.degree(),multiplicity) for p,multiplicity in mod2]==[(5,1)]
    assert sorted([(p.degree(),multiplicity) for p,multiplicity in mod5])==[(2,1),(3,1)]
    alpha=-a*v/4; beta=k*a*v**2/(8*b)
    span_det=sp.factor(sp.det(sp.Matrix([[alpha,-beta],[alpha,beta]])))
    assert sp.factor(span_det+k*a*a*v**3/(16*b))==0
    return {'result':'PASS','coefficients':'exact QQ[k,v,b,s,a,z]',
       'point_schemes':results,'generic_support':'line w=z=0 plus five points off the line',
       'separable_sample':{'k':1,'v':1,'b':'-1/4','s':'1/4','a':1,
           'P_plus':'z**5+z**2+1','mod_2_factors':[5], 'mod_5_factors':[2,3]},
       'last_span_determinant':str(span_det),
       'limitations':['Group-theoretic proof that the sample has no PGL2 symmetry is in paper',
                      'Generic non-periodicity is proved from point schemes, not from tangent-space rank',
                      'No assertion of distinct irreducible moduli components']}
if __name__=='__main__':
    print(json.dumps(verify(),indent=2,ensure_ascii=False))
