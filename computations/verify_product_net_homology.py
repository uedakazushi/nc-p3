#!/usr/bin/env python3
"""Generic eight-parameter product-net obstruction; dominance is checked separately."""
import argparse,json
from pathlib import Path
import sympy as s
import generic_boundary_field as field
from flint_rational_function import rational_function_field
from product_net_boundary import tensor,relation_windows,NAMES
from verify_product_net_dominance import dominance_certificate
from quintic_field_analysis import analyze_quintic,analyze_d2_normal_forms,analyze_d2_specialized_basis
from verify_nonkoszul_boundary import tensor_at


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--quotient-pivot',choices=('columns','full','full-batched'),default='full')
    parser.add_argument('--projection',choices=('guided','normal','dual'),default='guided')
    parser.add_argument('--polynomial-frame',action='store_true')
    parser.add_argument('--reference',choices=('generic','boolean'),default='generic')
    args=parser.parse_args()
    dominance=dominance_certificate()
    constants=rational_function_field(('z',))
    substitution=dict(zip(s.symbols(NAMES),(0,0,0,1,1,0,0,0)))
    class Specialized:
        zero,one=constants.zero,constants.one
        @staticmethod
        def from_sympy(expression):return constants.from_sympy(s.sympify(expression).subs(substitution))
    old=tensor(Specialized())
    expected=tensor_at(0)
    assert all(old.get(i,constants.zero).as_expr()==expected.get(i,0) for i in range(256))
    field.K,field.zero,field.one=constants,constants.zero,constants.one
    old_relations=relation_windows(old,Specialized())
    old_matrix=s.Matrix(16,6,lambda i,j:old_relations[1][j].get(divmod(i,4),constants.zero).as_expr())
    minor_rows=old_matrix.T.rref()[1]
    assert len(minor_rows)==6
    pair_minor=old_matrix[list(minor_rows),:].det()
    assert pair_minor!=0
    print('PASS dominance minor, old-point frame and generic pair-rank certificate',flush=True)
    K=rational_function_field(NAMES);field.K,field.zero,field.one=K,K.zero,K.one
    w=tensor(K,polynomial=args.polynomial_frame)
    print('PASS eight-variable cubic frame',flush=True)
    relations=relation_windows(w,K)
    assert len(relations[1])==6
    print('PASS generic pair rank6 from the nonzero old-point minor and the universal upper bound',flush=True)
    reference_parameters=(2,5,11,3,7,17,23,31) if args.reference=='generic' else (0,0,0,1,1,0,0,0)
    class PivotReference:
        zero,one=constants.zero,constants.one
        @staticmethod
        def from_sympy(expression):
            return constants.from_sympy(s.sympify(expression).subs(dict(zip(s.symbols(NAMES),reference_parameters))))
    pivot_relations=relation_windows(tensor(PivotReference()),PivotReference())
    if args.projection=='guided':report=analyze_d2_specialized_basis(relations,pivot_relations,K,constants)
    elif args.projection=='normal':report=analyze_d2_normal_forms(relations,K)
    else:report=analyze_quintic(relations,K,rank_strategy='full',d2_only=True,quotient_strategy=args.quotient_pivot)
    bound=120-report['d2_rank']-36
    certificate={'result':'PASS','arithmetic':'exact rational function field in eight independent variables',
        'parameters':list(NAMES),'polynomial_frame':args.polynomial_frame,'dominance_tangent_rank':dominance['tangent_rank'],
        'dominance_minor':dominance['nonzero_minor'],'old_point_frame_compatibility':True,
        'pair_rank_minor_rows':list(minor_rows),'pair_rank_minor_at_old_point':str(pair_minor),
        'pair_rank_reason':'old-point lower bound6 and universal annihilator upper bound6',
        'pivot_reference_parameters':list(reference_parameters),'projection':args.projection,'quotient_strategy':args.quotient_pivot,'quintic':report,'universal_d3_upper_bound':36,'generic_H2_lower_bound':bound,
        'obstructed_generic_chart':bound>0,
        'limits':['Transfer to every basepoint partial-symmetric candidate requires the universal-frame-bundle specialization proof in the manuscript.',
                  'No general low-rank AS classification follows from the product-net chart alone.']}
    if args.output:args.output.write_text(json.dumps(certificate,indent=2)+'\n')
    print(json.dumps(certificate,indent=2))


if __name__=='__main__':main()
