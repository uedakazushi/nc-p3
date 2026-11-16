#!/usr/bin/env python3
"""Exact first-order opposite-rank cut at the Clifford point; no higher hull claim."""
import argparse
import json
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from sympy import Matrix,zeros
from tensor_symmetrization import act_factor
from pbw_complex import clifford0,matrices,PAIRS,MON2,rank
from verify_moduli_jets import quotient_basis
from frobenius_tensor import nullspace_sparse,digits

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    rules=clifford0();N,M=matrices(rules)
    even=quotient_basis(N,M,1);odd=quotient_basis(N,M,-1)
    assert (even.cols,odd.cols)==(9,6)
    base_dual=[]
    for m in MON2:
        relation={m:F(1)}
        for a,b in PAIRS:
            coefficient=rules[b,a].get(m,0)
            if coefficient:relation[b,a]=F(coefficient)
        base_dual.append(relation)

    def idx(word):
        val=0
        for a in word:val=4*val+a
        return val
    rows=[];metadata=[]
    for position in range(3):
        for left in product(range(4),repeat=position):
            for right in product(range(4),repeat=2-position):
                for m,relation in enumerate(base_dual):
                    rows.append({idx(left+pair+right):c for pair,c in relation.items()})
                    metadata.append((position,left,right,m))
    w=nullspace_sparse(rows,256);assert len(w)==1
    w={i:c/w[0][27] for i,c in w[0].items()}
    rhs=[]
    for position,left,right,m in metadata:
        row=[]
        for j in range(15):
            vector=even[:,j] if j<9 else odd[:,j-9]
            scale=1 if j<9 else (-1)**position
            row.append(-sum(F(vector[10*r+m])*scale*w.get(idx(left+(b,a)+right),F(0)) for r,(a,b) in enumerate(PAIRS)))
        rhs.append(row)
    rows.append({27:F(1)});rhs.append([F(0)]*15)
    pivots={}
    for equation,values in zip(rows,rhs):
        eq=dict(equation);values=list(values)
        while eq:
            pivot=min(eq)
            if pivot not in pivots:
                scale=eq[pivot]
                pivots[pivot]=({i:c/scale for i,c in eq.items()},[c/scale for c in values]);break
            coefficients,known=pivots[pivot];scale=eq[pivot]
            for i,c in coefficients.items():
                value=eq.get(i,F(0))-scale*c
                if value:eq[i]=value
                else:eq.pop(i,None)
            values=[a-scale*b for a,b in zip(values,known)]
        else:
            assert not any(values),values
    assert len(pivots)==256
    solution={}
    for pivot,(eq,values) in sorted(pivots.items(),reverse=True):
        solution[pivot]=[value-sum(c*solution[i][j] for i,c in eq.items() if i!=pivot) for j,value in enumerate(values)]
    for equation,expected_values in zip(rows,rhs):
        assert all(sum(coefficient*solution[index][j] for index,coefficient in equation.items())==expected_values[j]
                   for j in range(15))
    base=Matrix(16,16,lambda row,col:w.get(idx((row//4,col//4,row%4,col%4)),F(0)))
    assert base.rank()==6 and base.T==base
    left=Matrix.hstack(*base.T.nullspace());right=Matrix.hstack(*base.nullspace())
    residual=[]
    for j in range(15):
        delta=Matrix(16,16,lambda row,col:solution[idx((row//4,col//4,row%4,col%4))][j])
        projected=left.T*delta*right
        residual.append(Matrix(100,1,list(projected)))
    condition=Matrix.hstack(*residual)
    assert condition[:,:9]==zeros(100,9) and rank(condition[:,9:])==6 and rank(condition)==6

    import sympy as s
    matrices_odd=[Matrix(10,10,list(condition[:,9+j])) for j in range(6)]
    assert all(a.T==-a for a in matrices_odd)
    stack=Matrix.vstack(*matrices_odd)
    active=stack.rref()[1]
    assert len(stack.nullspace())==6 and len(active)==4
    o=s.symbols('o0:6')
    A=sum((o[j]*matrices_odd[j].extract(active,active) for j in range(6)),zeros(len(active)))

    if len(active)==4:
        pf=s.expand(A[0,1]*A[2,3]-A[0,2]*A[1,3]+A[0,3]*A[1,2])
        Q=o[0]*o[5]+o[1]*o[2]-o[2]*o[3]+o[3]*o[4]
        assert s.expand(pf+Q)==0

    # Gauge directions do not contribute to the determinantal derivative.
    for factor in range(4):
        for column in range(16):
            elementary=zeros(4);elementary[column//4,column%4]=1
            variation=act_factor(w,elementary,factor)
            delta=Matrix(16,16,lambda row,col:variation.get(idx((row//4,col//4,row%4,col%4)),F(0)))
            assert left.T*delta*right==zeros(10)
    # The chosen even chart stays symmetric for all its formal parameters.
    assert all(even[10*r+m,j]==0 for r in range(6) for m,pair in enumerate(MON2)
               if pair[0]!=pair[1] for j in range(9))
    mode_ranks=[]
    for z in (s.Integer(1),s.Integer(-1),s.I,-s.I):
        r1=rank(N[0]+z*N[1]);r2=rank(M[0]+z*M[1])
        mode_ranks.append([str(z),r1,r2,60-r1-r2])
    assert mode_ranks==[['1',15,36,9],['-1',15,39,6],['I',15,45,0],['-I',15,45,0]]
    pivot_rows=condition[:,9:].T.rref()[1]
    minor=condition[:,9:].extract(pivot_rows,range(6))
    assert minor.det()!=0
    certificate={'result':'PASS','arithmetic':'exact rational and QQ(i)',
        'mode_ranks':mode_ranks,'flat_representation_tangent_dimension':75,
        'rank_six_representation_tangent_dimension':69,
        'based_rank_six_quotient_tangent_dimension':9,
        'even_derivative_zero':True,'odd_derivative_rank':6,
        'odd_derivative_minor_rows':list(pivot_rows),'odd_derivative_minor_det':str(minor.det()),
        'gauge_directions_checked':64,'top_tensor_lift_equations':481,
        'common_kernel_dimension':6,'active_indices':list(active),
        'active_skew_matrix':[[str(x) for x in A.row(i)] for i in range(4)],
        'pfaffian':str(pf),'comparison_with_known_obstruction':'Pfaffian=-Q',
        'limits':['Tangent/trace computation at the specified Clifford point; the finite-type local coverage uses the manuscript proof.',
                  'No global low-rank coverage or full Clifford deformation-ring calculation.',
                  'First-order low-rank constraints do not establish vanishing of all higher obstructions of the unrestricted problem.']}
    if args.output:
        args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
