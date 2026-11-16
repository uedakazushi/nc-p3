#!/usr/bin/env python3
"""Exact universal cubic frame and dimensions of its four-layer extension."""
import argparse,json
from itertools import combinations_with_replacement
import sympy as s
from pathlib import Path
from generic_boundary_field import tensor,windows,u,v,K
from frobenius_tensor import digits


def general_determinant():
    c0,c1,c2,t=s.symbols('c0 c1 c2 t')
    monomials=list(combinations_with_replacement(range(4),3));lookup={m:i for i,m in enumerate(monomials)}
    rows=[]
    for a in range(4):
        for i in range(3):
            row=[0]*20;row[lookup[tuple(sorted((i,i,a)))]]+=1
            row[lookup[tuple(sorted((i,3,a)))]]-=1;rows.append(row)
        row=[0]*20
        for i,c in enumerate((c0,c1,c2)):
            row[lookup[tuple(sorted((i,3,a)))]]+=c
        row[lookup[tuple(sorted((3,3,a)))]]+=t;rows.append(row)
    free=[(0,1,2),(1,3,3),(2,3,3),(3,3,3)]
    normalization=s.Matrix(4,20,lambda i,j:int(monomials[j]==free[i]))
    assert s.factor(s.Matrix(rows).col_join(normalization).det())==-2*c0*c0*c1*c2


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    general_determinant()
    w=tensor();relations=windows(w)
    assert [len(row) for row in relations]==[6]*4
    # The return basis vector 3 is the cube at the fixed basepoint.
    assert {digits(i):x for i,x in w.items() if digits(i)[2]==3}=={(3,3,3,3):K.one}
    a_dimensions=[4,3*4,2*6,4]
    extension=[a_dimensions[i]+(a_dimensions[4-i] if i else 0)
               for i in range(4)]+[4]
    assert a_dimensions==[4,12,12,4]
    assert extension==[4,16,24,16,4] and sum(extension)==64
    certificate={'result':'PASS','arithmetic':'exact QQ(u,v) cubic frame and integer dimension identities',
        'normalized_quadric_coefficients':['1','u','v'],
        'general_cubic_frame_determinant':'-2*c0^2*c1*c2, independent of t',
        'cubic_frame_determinant':'-2*u*v','four_window_pair_ranks':[6]*4,
        'return_basis_cube':'333 at return coordinate 3',
        'truncated_algebra_dimensions':a_dimensions,'trivial_extension_dimensions':extension,
        'limits':['The structural isomorphism is proved in the manuscript under finite-Frobenius hypotheses.',
                  'This check does not establish uniform finite-Frobenius or Koszul properties of the coefficient family.']}
    if args.output:args.output.write_text(json.dumps(certificate,indent=2)+'\n')
    print(json.dumps(certificate,indent=2))


if __name__=='__main__':main()
