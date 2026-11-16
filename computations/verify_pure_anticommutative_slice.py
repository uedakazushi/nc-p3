#!/usr/bin/env python3
"""Direct three-flattening tangent and CI-chart coverage at the pure anti point."""
import argparse
import json
from pathlib import Path
from fractions import Fraction as F
from itertools import product,combinations
from sympy import Matrix
from frobenius_tensor import top_tensor,digits,nullspace_sparse,independent_columns
from tensor_symmetrization import act_factor
from pbw_complex import PAIRS,clifford0


def encode(word):
    value=0
    for a in word:value=4*value+a
    return value


def from_rules(rules):
    relations=[]
    for a,b in PAIRS:
        relation={(b,a):F(1)}
        for pair,value in rules[b,a].items():relation[pair]=-F(value)
        relations.append(relation)
    w=top_tensor([relations]*4)
    scale=w[27]
    return {i:x/scale for i,x in w.items()}


def cut_equations(w):
    equations=[]
    for left in ((0,1),(1,2),(0,2)):
        right=tuple(i for i in range(4) if i not in left)
        matrix=Matrix.zeros(16)
        for i,value in w.items():
            word=digits(i)
            matrix[4*word[left[0]]+word[left[1]],4*word[right[0]]+word[right[1]]]=value
        assert matrix.rank()==6
        kl=Matrix.hstack(*matrix.T.nullspace());kr=Matrix.hstack(*matrix.nullspace())
        for a,b in product(range(10),repeat=2):
            row={}
            for i in range(256):
                word=digits(i)
                coefficient=F(kl[4*word[left[0]]+word[left[1]],a]
                              *kr[4*word[right[0]]+word[right[1]],b])
                if coefficient:row[i]=coefficient
            equations.append(row)
    return equations


def ci_variation(pair,square):
    # e_square^2 + t e_pair[0] e_pair[1] = 0 in the commutative dual.
    # At first order exactly one square can be replaced before taking trace.
    variation={}
    for word in product(range(4),repeat=4):
        if word.count(square)<2:continue
        remainder=list(word);remainder.remove(square);remainder.remove(square)
        if sorted(remainder+list(pair))==[0,1,2,3]:variation[encode(word)]=F(-1)
    return variation


def projective(tensor,w):
    scale=tensor.get(27,F(0))
    result={i:tensor.get(i,F(0))-scale*w.get(i,F(0)) for i in tensor.keys()|w.keys()}
    return {i:x for i,x in result.items() if x}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    pure={(b,a):{(a,b):-1} for a,b in PAIRS}
    w=from_rules(pure)
    assert w=={encode(word):F(1) for word in product(range(4),repeat=4) if len(set(word))==4}
    equations=cut_equations(w)
    assert len(equations)==300
    kernel=nullspace_sparse(equations,256)
    assert len(kernel)==70
    orbit=[]
    for factor in range(4):
        for index in range(16):
            matrix=Matrix.zeros(4);matrix[index//4,index%4]=1
            orbit.append(projective(act_factor(w,matrix,factor),w))
    ci=[ci_variation(pair,square) for square in range(4) for pair in PAIRS]
    for vector in orbit+ci:
        assert all(sum(value*vector.get(i,F(0)) for i,value in row.items())==0 for row in equations)
    orbit_rank=len(independent_columns(orbit));ci_rank=len(independent_columns(ci))
    combined=len(independent_columns(orbit+ci))
    assert (orbit_rank,ci_rank,combined)==(57,24,69)
    # Independent cross-check at the previously studied Clifford point.
    other=from_rules(clifford0())
    other_kernel=nullspace_sparse(cut_equations(other),256)
    assert len(other_kernel)==70
    # Alternating determinant tensor: the product-group orbit already
    # fills the whole direct three-cut tangent space.
    alternating=from_rules({(b,a):{(a,b):1} for a,b in PAIRS})
    alternating_equations=cut_equations(alternating)
    alternating_kernel=nullspace_sparse(alternating_equations,256)
    assert len(alternating_kernel)==46
    alternating_orbit=[]
    for factor in range(4):
        for index in range(16):
            matrix=Matrix.zeros(4);matrix[index//4,index%4]=1
            vector=projective(act_factor(alternating,matrix,factor),alternating)
            assert all(sum(value*vector.get(i,F(0)) for i,value in row.items())==0
                       for row in alternating_equations)
            alternating_orbit.append(vector)
    assert len(independent_columns(alternating_orbit))==45
    certificate={'result':'PASS','arithmetic':'exact rational',
        'three_cut_equations':300,'affine_tangent_dimension':70,
        'projective_tangent_dimension':69,'projective_orbit_tangent_rank':57,
        'symmetric_CI_chart_tangent_rank':24,'combined_smooth_source_rank':69,
        'quotient_first_order_dimension':12,
        'Clifford_C0_direct_tangent_crosscheck':69,
        'alternating_projective_tangent_dimension':45,
        'alternating_projective_orbit_rank':45,
        'limits':['69 is a presentation dimension, 12 is a first-order quotient at a stabilizer-jumping point.',
                  'No 12-dimensional component or global low-rank coverage is asserted.',
                  'Open-neighborhood coverage uses the smooth-source proof in the manuscript.']}
    if args.output:args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
