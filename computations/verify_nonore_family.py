#!/usr/bin/env python3
"""Exact all-parameter PBW and moduli checks for a non-Ore 2-period family.

Uses exact QQ[k,v,b,s,a] polynomials and rational tangent matrices.
All eight degree-three ambiguities are checked, with a well-founded
positive weight order for every rewrite. AS regularity is the manuscript's
Koszul-Frobenius filtered argument; no global period-one exclusion claimed.
"""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import argparse
import json
import sympy as sp
from pbw_complex import PAIRS, MON2, clifford0, matrices, rank

k,v,b,s,a=sp.symbols('k v b s a')
PARAMETERS=(k,v,b,s,a)
BASE={k:1,v:1,b:1,s:1,a:0}
WEIGHTS=(3,5,2,2)  # W,Y,Z,X


def rules(epsilon=1):
    e=sp.sympify(epsilon)
    return {
        (1,0):{(0,1):-1,(2,3):e*a},
        (2,0):{(0,2):-1},
        (2,1):{(1,2):-1,(0,0):k},
        (3,0):{(0,3):-1,(2,2):v},
        (3,1):{(1,3):-1,(0,0):b,(2,2):s},
        (3,2):{(2,3):-1},
    }


def overlaps():
    rr=(rules(1),rules(-1))
    def ordered(word):
        return (sum(WEIGHTS[j] for j in word),word)
    for block in rr:
        for lhs,rhs in block.items():
            assert all(ordered(word)<ordered(lhs) for word in rhs)
    @lru_cache(None)
    def normal(word,start):
        for pos in range(len(word)-1):
            if word[pos]>word[pos+1]:
                out={}
                for w,co in rr[(start+pos)%2][word[pos:pos+2]].items():
                    for result, value in normal(word[:pos]+w+word[pos+2:],start).items():
                        out[result]=sp.expand(out.get(result,0)+co*value)
                return {w:c for w,c in out.items() if c!=0}
        return {word:sp.Integer(1)}
    total=0
    for start in (0,1):
        for tri in combinations(range(4),3):
            word=tuple(reversed(tri))
            difference={}
            for pos,sign in ((0,1),(1,-1)):
                for repl, co in rr[(start+pos)%2][word[pos:pos+2]].items():
                    for out, val in normal(word[:pos]+repl+word[pos+2:],start).items():
                        difference[out]=sp.expand(difference.get(out,0)+sign*co*val)
            assert all(sp.expand(x)==0 for x in difference.values()),(start,word,difference)
            total+=1
    assert total==8
    return total


def tangent_and_torus():
    N,M=matrices(clifford0())
    rs=rules(1)
    A=sp.Matrix([rs[(hi,lo)].get(mon,0) for lo,hi in PAIRS for mon in MON2])
    J=A.jacobian(PARAMETERS).subs(BASE)
    e=J[:,:4]; o=J[:,4:]
    assert (M[0]+M[1])*e==sp.zeros(80,4)
    assert (M[0]-M[1])*o==sp.zeros(80,1)
    r_orbit=rank(N[0].row_join(N[1]))
    r_params=rank(N[0].row_join(N[1]).row_join(J))
    assert (r_orbit,r_params)==(30,32)
    assert rank((N[0]+N[1]).row_join(e))==16
    assert rank((N[0]-N[1]).row_join(o))==16
    # Constant diagonal scaling weights for (k,v,b,s,a).
    scaling=sp.Matrix([[2,-1,-1,0],[-1,0,2,-1],
                       [2,-1,0,-1],[0,-1,2,-1],[-1,-1,1,1]])
    assert scaling.rank()==3
    inv1=sp.Matrix([-2,-2,1,1,0])
    inv2=sp.Matrix([-5,-3,4,0,1])
    assert inv1.T*scaling==sp.zeros(1,4)
    assert inv2.T*scaling==sp.zeros(1,4)
    return {'single_relation_orbit_rank':r_orbit,
            'orbit_plus_family_rank':r_params,
            'even_tangent_class_dimension':1,
            'odd_tangent_class_dimension':1,
            'torus_effective_rank':3,
            'rational_scaling_invariants':['b*s/(k**2*v**2)', 'a*b**4/(k**5*v**3)'],
            'indexed_moduli_image_dimension':2}


def verify():
    return {'result':'PASS',
            'base':'A5_k,v,b,s,a (all geometric fibers)',
            'arithmetic':'exact QQ[k,v,b,s,a] and rational matrices',
            'weights_W_Y_Z_X':list(WEIGHTS),
            'polynomial_overlaps':overlaps(),
            'moduli':tangent_and_torus(),
            'not_claimed':['generic minimal period two',
                'full AS regular Z-algebra component',
                'that all fibers are outside the existing nine-parameter image',
                'formal or peer review']}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=verify()
    text=json.dumps(result,indent=2,ensure_ascii=False)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    print(text)
