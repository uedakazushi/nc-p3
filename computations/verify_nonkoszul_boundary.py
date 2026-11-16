#!/usr/bin/env python3
"""Finite-Frobenius boundary of the symmetric chart with first defect in degree five."""
import argparse
import json
from pathlib import Path
from itertools import combinations_with_replacement
from fractions import Fraction as F
import sympy as s
from frobenius_tensor import (digits,flatten,reconstructed_windows,orthogonal_windows,
                             insertion_columns,nullspace_sparse,exact_rank)
from tensor_symmetrization import symmetrizer_space,symmetrize
from verify_quartic_obstruction import quintic_d3_rank


def cubic_frame():
    t=s.symbols('t');monomials=list(combinations_with_replacement(range(4),3))
    lookup={m:i for i,m in enumerate(monomials)};rows=[]
    for a in range(4):
        for i in range(3):
            row=[0]*20
            row[lookup[tuple(sorted((i,i,a)))]]+=1
            row[lookup[tuple(sorted((i,3,a)))]]-=1
            rows.append(row)
        row=[0]*20
        for i in range(3):row[lookup[tuple(sorted((i,3,a)))]]+=1
        row[lookup[tuple(sorted((3,3,a)))]]+=t
        rows.append(row)
    equations=s.Matrix(rows)
    free=[(0,1,2),(1,3,3),(2,3,3),(3,3,3)]
    normalization=s.Matrix(4,20,lambda i,j:int(monomials[j]==free[i]))
    square=equations.col_join(normalization)
    assert s.expand(square.det())==-2
    basis=square.inv()*s.zeros(16,4).col_join(s.eye(4))
    assert equations*basis==s.zeros(16,4)
    assert normalization*basis==s.eye(4)
    return t,monomials,basis


def tensor_at(parameter):
    t,monomials,basis=cubic_frame();lookup={m:i for i,m in enumerate(monomials)}
    w={}
    for encoded in range(256):
        a,b,j,c=digits(encoded)
        coefficient=F(basis[lookup[tuple(sorted((a,b,c)))],j].subs(t,parameter))
        if coefficient:w[encoded]=coefficient
    return w


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    w=tensor_at(0)
    assert len(w)==83
    assert [flatten(w,(i,)) for i in range(4)]==[4]*4
    assert [flatten(w,(0,j)) for j in (1,2,3)]==[6]*3
    windows=reconstructed_windows(w);dual=orthogonal_windows(windows)
    intersections=[[len(nullspace_sparse(insertion_columns(dual,i,n),4**n)) for n in (3,4)] for i in range(4)]
    assert intersections==[[4,1]]*4
    assert all(len(nullspace_sparse(insertion_columns(dual,i,5),1024))==0 for i in range(4))
    spaces=[symmetrizer_space(w,j) for j in (1,2,3)]
    assert spaces[0]==spaces[2]==[s.eye(4)]
    assert spaces[1]==[s.diag(0,0,0,1)]
    assert symmetrize(w)[0] is None
    assert symmetrize(w,-1)[0] is None
    quartic=[256-exact_rank(insertion_columns(windows,i,4)) for i in range(4)]
    quintic=[1024-exact_rank(insertion_columns(windows,i,5)) for i in range(4)]
    assert quartic==[35]*4 and quintic==[63,56,56,56]
    d3_ranks=[quintic_d3_rank(windows,i) for i in range(4)]
    assert d3_ranks==[36]*4
    generic=tensor_at(1)
    assert symmetrize(generic)[0] is not None
    certificate={'result':'PASS','arithmetic':'exact rational polynomial and linear algebra',
        'quadrics':['e0^2-e0*e3','e1^2-e1*e3','e2^2-e2*e3','e3*(e0+e1+e2+t*e3)'],
        'cubic_normalization_coordinates':['012','133','233','333'],
        'cubic_normalization_determinant':-2,
        'symmetric_CI_parameters':'t not in {0,-1,-2,-3}',
        'special_tensor_terms':83,'singleton_ranks':[4]*4,'pair_cut_ranks':[6]*3,
        'intersection_dimensions_W3_W4':intersections,'dual_degree_five_dimensions':[0]*4,
        'symmetrizers':['k I','k E33','k I'],'alternating_membership':False,
        'primal_degree_four_dimensions':quartic,'primal_degree_five_dimensions':quintic,
        'quintic_d3_ranks':d3_ranks,'quintic_H2_dimensions':[7,0,0,0],
        'quintic_H3_dimensions':[0]*4,
        'primal_relation_insertion_rank_degree_five':[1024-d for d in quintic],
        'conclusions':['The symmetric chart is not relatively closed inside the finite-Frobenius low-rank candidate locus.',
                       'Finite Frobenius plus all four quartic Koszul filters does not imply Koszulity.',
                       'The special point is not a standard AS object and is not a LOW-UNCLASSIFIED AS witness.'],
        'limits':['No full AS classification or exhaustive analysis of the symmetric boundary.','Generic CI fibers are regular by the existing theorem; special-fiber non-Koszulity is certified by a degree-five Hilbert defect.']}
    if args.output:args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
