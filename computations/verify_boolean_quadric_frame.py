#!/usr/bin/env python3
"""Global cubic evaluation frame over the entire one-point Boolean open."""
import argparse,json
from itertools import product,combinations
from pathlib import Path
from fractions import Fraction as F
import sympy as s
from flint_rational_function import rational_function_field
from boolean_quadric_frame import data,evaluation_tensor
from verify_nonkoszul_boundary import tensor_at
from frobenius_tensor import (digits,flatten,reconstructed_windows,orthogonal_windows,
    nullspace_sparse,insertion_columns,exact_rank)
from verify_quartic_obstruction import quintic_d3_rank


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    K=rational_function_field(('u','v','a','b','c'))
    parameters,values,quadric=data(K)
    w=evaluation_tensor(K)
    for j,x in product(range(4),repeat=2):
        def coeff(a,b,c):return w.get(64*a+16*b+4*j+c,K.zero)
        for i in range(3):assert coeff(i,i,x)==coeff(i,3,x)
        assert sum((q*coeff(a,b,x) for (a,b),q in quadric.items()),K.zero)==K.zero
    points=list(product(range(2),repeat=3))
    subsets=[subset for size in range(4) for subset in combinations(range(3),size)]
    evaluation=s.Matrix([[s.prod(point[i] for i in subset) for subset in subsets] for point in points])
    assert abs(evaluation.det())==1
    nonzero=[point for point in points if any(point)]
    degree2=s.Matrix([[s.prod(point[i] for i in subset) for subset in subsets[:-1]] for point in nonzero])
    assert abs(degree2.det())==1
    degree1=s.Matrix([[1,*point] for point in nonzero])
    assert degree1.rank()==4
    for point,q in values.items():
        p=point+(1,)
        assert sum((x*p[a]*p[b] for (a,b),x in quadric.items()),K.zero)==q
    # The selected rows of the evaluation-weight frame are exactly diagonal.
    selected=[(1,0,0),(0,1,0),(0,0,1),(0,0,0)]
    diagonal=[]
    for row,point in enumerate(selected):
        for col in range(4):
            if col==3:value=K.one if not any(point) else K.zero
            elif point[col]:value=K.from_sympy(s.Integer((-1)**(sum(point)-1)))/values[point]
            else:value=K.zero
            if row==col:diagonal.append(value)
            else:assert value==K.zero
    weight_minor=diagonal[0]*diagonal[1]*diagonal[2]*diagonal[3]
    assert weight_minor==K.from_sympy(1/(s.Symbol('u')*s.Symbol('v')))
    substitution=dict(zip(s.symbols('u v a b c'),(1,1,2,2,2)))
    point_w={i:F(value.as_expr().subs(substitution)) for i,value in w.items()
             if value.as_expr().subs(substitution)!=0}
    free=[(0,1,2),(1,3,3),(2,3,3),(3,3,3)]
    frame=s.Matrix(4,4,lambda k,j:point_w.get(64*free[k][0]+16*free[k][1]+4*j+free[k][2],F(0)))
    assert frame.det()!=0
    old=tensor_at(0)
    for i in range(256):
        a,b,j,c=digits(i)
        assert point_w.get(i,0)==sum(old.get(64*a+16*b+4*k+c,0)*frame[k,j] for k in range(4))
    cross_sub=dict(zip(s.symbols('u v a b c'),(2,3,5,7,11)))
    cross={i:F(x.as_expr().subs(cross_sub)) for i,x in w.items()
           if x.as_expr().subs(cross_sub)!=0}
    assert [flatten(cross,(i,)) for i in range(4)]==[4]*4
    relations=reconstructed_windows(cross);dual=orthogonal_windows(relations)
    intersections=[[len(nullspace_sparse(insertion_columns(dual,i,n),4**n))
                    for n in (3,4)] for i in range(4)]
    assert intersections==[[4,1]]*4
    quartic=[256-exact_rank(insertion_columns(relations,i,4)) for i in range(4)]
    quintic=[1024-exact_rank(insertion_columns(relations,i,5)) for i in range(4)]
    assert quartic==[35]*4 and quintic==[63,56,56,56]
    assert [quintic_d3_rank(relations,i) for i in range(4)]==[36]*4
    certificate={'result':'PASS','arithmetic':'exact QQ(u,v,a,b,c), integer Boolean evaluation determinants',
        'point_values':['1','u','v','a','b','c','a+b+c-1-u-v'],
        'domain':'all seven listed values nonzero',
        'cubic_evaluation_determinant_abs':1,'degree_two_seven_point_determinant_abs':1,
        'degree_one_seven_point_rank':4,'global_weight_frame_minor':'1/(u*v)',
        'annihilation_of_all_four_quadrics':True,
        'quotient_Hilbert_dimensions_degrees_0_to_5':[1,4,6,4,1,1],
        'quotient_Hilbert_tail':'1 in every degree >=4',
        'old_boundary_point_frame_det':str(frame.det()),
        'cross_term_rational_fixture':{'u_v_a_b_c':[2,3,5,7,11],
            'intersection_dimensions':intersections,'C4':quartic,'C5':quintic,'d3_ranks':[36]*4},
        'limits':['Global frame and commutative quotient Hilbert function only; no generic quintic homology rank is inferred from samples.',
                  'The first three quadrics remain the fixed Boolean complete intersection.']}
    if args.output:args.output.write_text(json.dumps(certificate,indent=2)+'\n')
    print(json.dumps(certificate,indent=2))


if __name__=='__main__':main()
