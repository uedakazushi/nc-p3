#!/usr/bin/env python3
"""Exact non-Frobenius boundary: ranks stay six but symmetrizers become singular."""
import argparse
import json
from pathlib import Path
from fractions import Fraction as F
from sympy import Matrix,symbols
from verify_ci_centroid import from_rules
from pbw_complex import PAIRS
from frobenius_tensor import (digits,flatten,reconstructed_windows,orthogonal_windows,
                             nullspace_sparse,insertion_columns)
from tensor_symmetrization import symmetrizer_space,invertible_witness


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    rules={(b,a):{(a,b):-1} for a,b in PAIRS};rules[1,0][3,3]=1
    source=from_rules(rules)
    weights=[[3,3,1,0],[3,3,1,0],[0,0,4,3],[3,3,1,0]]
    costs={i:sum(weights[j][a] for j,a in enumerate(digits(i))) for i in source}
    assert min(costs.values())==4
    boundary={i:x for i,x in source.items() if costs[i]==4}
    expected={}
    for i in range(256):
        word=digits(i)
        if len(set(word))==4 and word[2] in (0,1):expected[i]=F(1)
        elif sorted(word)==[2,3,3,3]:expected[i]=F(-1)
    assert boundary==expected
    assert all(flatten(boundary,(j,))==4 for j in range(4))
    assert all(flatten(boundary,(0,j))==6 for j in (1,2,3))
    spaces=[symmetrizer_space(boundary,j) for j in (1,2,3)]
    assert spaces[0]==spaces[2]==[Matrix.eye(4)]
    assert len(spaces[1])==2
    e32=Matrix.zeros(4);e32[3,2]=1
    e22e33=Matrix.diag(0,0,1,1)
    assert spaces[1]==[e32,e22e33]
    a,b=symbols('a b')
    assert (a*e32+b*e22e33).det()==0
    assert invertible_witness(spaces[1]) is None
    windows=reconstructed_windows(boundary);dual=orthogonal_windows(windows)
    dimensions=[[len(nullspace_sparse(insertion_columns(dual,i,n),4**n)) for n in (3,4)] for i in range(4)]
    assert dimensions==[[6,10],[10,10],[6,6],[6,6]]
    basis=nullspace_sparse(insertion_columns(dual,0,3),64)
    def encode(word):return sum(int(c)*4**(2-j) for j,c in enumerate(word))
    expected_basis=[{'320':1,'230':1},{'321':1,'231':1},{'330':1},{'331':1},
                    {'332':1,'323':1,'233':1,'210':-1,'201':-1,'120':-1,'021':-1},
                    {'333':1,'310':-1,'301':-1,'130':-1,'031':-1}]
    assert basis==[{encode(word):F(value) for word,value in row.items()} for row in expected_basis]
    certificate={'result':'PASS','arithmetic':'exact rational linear algebra and integer weights',
                 'source':'D1 symmetric CI tensor from the regular specialization family',
                 'factor_weights':weights,'minimum_weight':4,'boundary_terms':16,
                 'singleton_ranks':[4]*4,'pair_cut_ranks':[6]*3,
                 'symmetrizer_dimensions':[1,2,1],
                 'middle_symmetrizers':'span(E32, E22+E33), determinant identically zero',
                 'intersection_dimensions_W3_W4':dimensions,'W03_basis':expected_basis,
                 'limits':['The point lies in Z_sym but fails the finite-Frobenius locus.',
                           'It does not answer whether the symmetric image is relatively closed within the finite-Frobenius locus.',
                           'All flattening ranks and all nonzero symmetrizer spaces alone are insufficient to certify Frobenius or AS regularity.']}
    if args.output:args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
