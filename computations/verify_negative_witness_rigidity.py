#!/usr/bin/env python3
"""Exact checks for partial-alternation normal forms and the Clifford model."""
import argparse,json
from pathlib import Path
from itertools import combinations
import sympy as s
from frobenius_tensor import digits
from negative_witness_normal_forms import tensor,pfaffian_matrices
from flint_rational_function import rational_function_field
import generic_boundary_field as field


def cut_rank(w,left):
    right=[j for j in range(4) if j not in left];matrix=s.zeros(16)
    for encoded,value in w.items():
        word=digits(encoded);matrix[4*word[left[0]]+word[left[1]],4*word[right[0]]+word[right[1]]]=value
    return matrix.rank()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
    variables,X,Y,q=pfaffian_matrices();assert (X*Y).applyfunc(s.expand)==q*s.eye(4)
    assert (Y*X).applyfunc(s.expand)==q*s.eye(4)
    Gamma=s.zeros(8);Gamma[:4,4:]=X;Gamma[4:,:4]=Y
    generators=[Gamma.diff(z) for z in variables]
    products=[]
    for mask in range(64):
        product=s.eye(8)
        for j,g in enumerate(generators):
            if mask>>j&1:product=product*g
        products.append(list(product))
    clifford_rank=s.Matrix.hstack(*(s.Matrix(p) for p in products)).rank();assert clifford_rank==64
    X_tensor={}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    value=int(a==c and b==d)-int(a==d and b==c)
                    if value:X_tensor[64*a+16*b+4*c+d]=value
    assert cut_rank(X_tensor,(0,2))==16
    K=rational_function_field(('z',));field.K,field.zero,field.one=K,K.zero,K.one
    reports={}
    for kind in ('volume','rank4','rank5'):
        w=tensor(kind);ranks=[cut_rank(w,pair) for pair in [(0,1),(1,2),(0,2)]]
        relations=field.windows({i:K.from_sympy(x) for i,x in w.items()});dual=field.dual_windows(relations)
        intersections=[(len(field.nullspace(field.inserted(dual,i,3),64)),
                        len(field.nullspace(field.inserted(dual,i,4),256))) for i in range(4)]
        singleton=[]
        for i in range(4):
            mat=s.zeros(4,64)
            for encoded,value in w.items():
                word=digits(encoded);rest=[word[j] for j in range(4) if j!=i]
                mat[word[i],16*rest[0]+4*rest[1]+rest[2]]=value
            singleton.append(mat.rank())
        assert singleton==[4]*4
        if kind=='volume':assert ranks==[6,6,6] and intersections==[(4,1)]*4
        else:assert ranks==[6,6,int(kind[-1])] and intersections==[(10,16),(8,16),(10,16),(8,16)]
        basis=field.nullspace(field.inserted(dual,0,3),64)
        reports[kind]={'singleton_ranks':singleton,'adjacent_then_opposite_ranks':ranks,
                       'W3_W4':intersections,'W3_start0_basis':[{''.join(map(str,digits(i,3))):str(x.as_expr()) for i,x in v.items()} for v in basis],
                       'tensor':{''.join(map(str,digits(i))):x for i,x in sorted(w.items())}}
    result={'result':'PASS','arithmetic':'exact QQ and polynomial identities',
            'pfaffian_factorization_identity':True,'Clifford_generated_matrix_algebra_dimension':clifford_rank,
            'identity_wedge_opposite_rank':16,'normal_forms':reports,
            'limits':['Normal-form exhaustiveness is proved in the manuscript; it is not inferred from these samples.','No claim of full AS classification.']}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
