#!/usr/bin/env python3
"""Exact quartic one-return Tor fixtures: CI, two points, and a tangent jet."""
import argparse,json
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from flint_rational_function import rational_function_field
from frobenius_tensor import digits
from trivial_extension_bitor import one_return_complex
from quintic_field_analysis import combine
import generic_boundary_field as field


def analyze(quadrics):
    K=rational_function_field(('z',));field.K,field.zero,field.one=K,K.zero,K.one
    pairs=list(combinations_with_replacement(range(4),2));triples=list(combinations_with_replacement(range(4),3))
    matrix=s.Matrix([[q.get(pair,0) for pair in pairs] for q in quadrics]);assert matrix.rank()==4
    free2=[pair for j,pair in enumerate(pairs) if j not in matrix.rref()[1]]
    b2=[]
    for v in matrix.nullspace():
        row={}
        for pair,value in zip(pairs,v):
            if value:row[pair]=row[pair[::-1]]=K.from_sympy(value)
        b2.append(row)
    cubic=[]
    for q in quadrics:
        for i in range(4):
            row=[0]*20
            for pair,value in q.items():row[triples.index(tuple(sorted(pair+(i,))))]+=value
            cubic.append(row)
    b3=s.Matrix(cubic).nullspace();assert len(b3)==4
    w={}
    for n in range(256):
        a,b,j,c=digits(n);value=b3[j][triples.index(tuple(sorted((a,b,c))))]
        if value:w[n]=K.from_sympy(value)
    exact_quadrics=[{pair:K.from_sympy(x) for pair,x in q.items()} for q in quadrics]
    d,bases=one_return_complex(w,b2,exact_quadrics,K,coordinate_pairs=free2,start_position=1,end_position=1)
    assert all(not combine(d[1],v,K.zero) for v in d[2])
    assert all(not combine(d[2],v,K.zero) for v in d[3])
    ranks=[field.sparse_full_pivot_rank(d[h]) for h in (1,2,3)]
    homology=[1-ranks[0],32-ranks[0]-ranks[1],156-ranks[1]-ranks[2],160-ranks[2]]
    return {'complex_dimensions':[1,32,156,160],'ranks':ranks,'homology':homology,'direct_zero_compositions':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
    boolean=[{(i,i):1,(i,3):-1} for i in range(3)]
    ci=boolean+[{(i,3):1 for i in range(4)}]
    two=boolean+[{(1,3):1,(2,3):2,(0,1):1,(0,2):1}]
    fat=[{(i,i):1} for i in range(3)]+[{(0,3):1,(1,2):1}]
    reports={name:analyze(q) for name,q in [('complete_intersection',ci),('two_distinct_points',two),('nonreduced_point',fat)]}
    assert reports['complete_intersection']['homology']==[0,0,0,35]
    assert reports['two_distinct_points']['homology'][3]>=42
    assert reports['nonreduced_point']['homology'][3]>=42
    # Verify the displayed length-two subschemes directly.
    epsilon=s.symbols('epsilon')
    for q in two:
        for p in [(0,0,0,1),(1,0,0,1)]:assert sum(x*p[a]*p[b] for (a,b),x in q.items())==0
    for q in fat:
        p=(0,epsilon,0,1);assert s.rem(sum(x*p[a]*p[b] for (a,b),x in q.items()),epsilon**2,epsilon)==0
    result={'result':'PASS','arithmetic':'exact QQ','fixtures':reports,
            'two_point_quadric':'x1*h+2*x2*h+x0*x1+x0*x2',
            'nonreduced_quadrics':['x0^2','x1^2','x2^2','x0*h+x1*x2'],
            'length_two_witnesses_checked':True,
            'limits':['The universal result is proved using the point/jet subcomplex, not inferred from these fixtures.','No finite-Frobenius property is asserted for the two obstructed fixtures.']}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
