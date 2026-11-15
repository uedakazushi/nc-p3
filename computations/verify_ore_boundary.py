#!/usr/bin/env python3
"""Exact PBW boundary and singular orbit-closure certificate (no numerics)."""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import argparse, json
import sympy as sp
from pbw_complex import PAIRS, MON2, add, clifford0, matrices, rank
from verify_algebraic_family import PARAMETERS, BASE_POINT, rules as seven

K,V,B,S,A,D=sp.symbols('K V B S A D')
weights=(1,2,2,4)  # primed generators W,Y,Z,T

def rules(e):
    return {(1,0):{(0,1):-1}, (2,0):{(0,2):-1},
      (2,1):{(1,2):-1,(0,0):1},
      (3,0):{(0,3):-1,(1,1):K},
      (3,1):{(1,3):-1,(0,0):D,(2,2):D,(0,2):e*A*V},
      (3,2):{(2,3):-1,(1,1):B/V,(0,0):S/V}}

def reduce_word(rs, word, start=0):
    @lru_cache(None)
    def nf(w,i):
        for j in range(len(w)-1):
            if w[j]>w[j+1]:
                out={}
                for mon,c in rs[(i+j)%2][w[j:j+2]].items():
                    for term,v in nf(w[:j]+mon+w[j+2:],i).items():
                        add(out,term,sp.expand(c*v))
                return {m:sp.factor(c) for m,c in out.items() if sp.factor(c)!=0}
        return {w:sp.Integer(1)}
    return nf(word,start)

def reduce_poly(rs, poly, start=0):
    out={}
    for word,c in poly.items():
        for term,v in reduce_word(rs,word,start).items():
            add(out,term,sp.expand(c*v))
    return {term:sp.factor(c) for term,c in out.items() if sp.factor(c)!=0}

def verify_pbw():
    rs=[rules(1),rules(-1)];n=0
    for start in (0,1):
        for pair,rhs in rs[start].items():
            assert all((sum(weights[j] for j in w),w)<
                 (sum(weights[j] for j in pair),pair) for w in rhs)
        for triple in combinations(range(4),3):
            word=tuple(reversed(triple));diff={}
            for j,sign in ((0,1),(1,-1)):
                for mon,c in rs[(start+j)%2][word[j:j+2]].items():
                    for term,v in reduce_word(rs,word[:j]+mon+word[j+2:],start).items():
                        add(diff,term,sp.expand(sign*c*v))
            assert all(sp.factor(c)==0 for c in diff.values()),(start,word,diff)
            n+=1
    return n

def verify_isomorphism():
    # Original generators W,Y,Z,X = (Y_primed,T_primed,W_primed,V Z_primed).
    rs=[rules(1),rules(-1)];source={0:(1,1),1:(3,1),2:(0,1),3:(2,V)}
    original=[
      {(1,0):1,(0,1):1,(2,3):-A,(2,2):-D,(3,3):-D/V**2},
      {(2,0):1,(0,2):1},
      {(2,1):1,(1,2):1,(0,0):-K},
      {(3,0):1,(0,3):1,(2,2):-V},
      {(3,1):1,(1,3):1,(0,0):-B,(2,2):-S},
      {(3,2):1,(2,3):1}]
    for start in (0,1):
        for rr in original:
            poly={}
            for word,c in rr.items():
                (a,aa),(b,bb)=source[word[0]],source[word[1]]
                coefficient=c*(-1 if start and word==(2,3) and rr is original[0] else 1)
                add(poly,(a,b),coefficient*aa*bb)
            assert reduce_poly(rs,poly,start)=={}
        # For D != 0 use Tnew = T/D - S*Y/(V*D).
        # The old seven-family coefficients are (K/D,0,-2S/(VD),B/(VD),0,AV/D,0).
        vals=dict(zip(PARAMETERS,(K/D,0,-2*S/(V*D),B/(V*D),0,A*V/D,0)))
        Tnew={(3,):1/D,(1,):-S/(V*D)}
        for b in (0,1,2):
            poly={}
            for a,ca in Tnew.items():
                add(poly,a+(b,),ca);add(poly,(b,)+a,ca)
            target=seven(epsilon=1 if start==0 else -1)[(3,b)]
            for word,c in target.items():
                if word != (b,3): add(poly,word,-sp.sympify(c).subs(vals))
            assert reduce_poly(rs,poly,start)=={}
    return {'seven_parameters':{str(k):str(v) for k,v in vals.items()},
            'change_of_generators':['W_new=Z','Y_new=W','Z_new=X/V',
                'X_new=Y/D - S*W/(V*D)']}

def verify_tangent():
    N,_=matrices(clifford0());rr=seven()
    coeff=sp.Matrix([rr.get((hi,lo),{}).get(mon,0) for lo,hi in PAIRS for mon in MON2])
    J=coeff.jacobian(PARAMETERS).subs(dict(zip(PARAMETERS,BASE_POINT)))
    d=sp.Matrix([int(pair==(0,1) and mon in ((2,2),(3,3))) for pair in PAIRS for mon in MON2])
    o=sp.Matrix([int(pair==(0,1) and mon==(2,3)) for pair in PAIRS for mon in MON2])
    orbit=N[0].row_join(N[1]);base=orbit.row_join(J)
    ranks=[rank(orbit),rank(base),rank(base.row_join(d)),rank(base.row_join(o)),
           rank(base.row_join(d).row_join(o))]
    assert ranks==[30,37,38,38,39],ranks
    return {'orbit_and_boundary_ranks':ranks,'orbit_closure_dimension':37,
            'tangent_lower_bound_at_Clifford':39}

def verify():
    return {'result':'PASS','exact_arithmetic':'QQ(K,V,B,S,A,D)',
       'parameter_base':'Spec k[K,V,V^-1,B,S,A,D]',
       'two_periodic_PBW_overlaps':verify_pbw(),
       'all_six_original_relations_and_generic_seven_isomorphism':verify_isomorphism(),
       'singularity_certificate':verify_tangent(),
       'scope':['NonOre five-family image is in seven-family orbit-saturated closure',
                'That 37-dimensional Grassmannian image closure is singular at Clifford point',
                'No assertion about distinct full moduli irreducible components']}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path)
    a=p.parse_args();out=json.dumps(verify(),indent=2,ensure_ascii=False)+'\n'
    if a.output:a.output.write_text(out)
    print(out)
