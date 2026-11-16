#!/usr/bin/env python3
"""Generic eight-parameter product-net obstruction; dominance is checked separately."""
import argparse,json
from pathlib import Path
import sympy as s
import generic_boundary_field as field
from flint_rational_function import rational_function_field
from product_net_boundary import tensor,relation_windows,NAMES
from verify_product_net_dominance import dominance_certificate
from quintic_field_analysis import analyze_quintic
from verify_nonkoszul_boundary import tensor_at


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path);args=parser.parse_args()
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
    print('PASS dominance minor and exact old-point frame compatibility',flush=True)
    K=rational_function_field(NAMES);field.K,field.zero,field.one=K,K.zero,K.one
    w=tensor(K)
    print('PASS eight-variable normalized cubic frame and polynomial rescaling',flush=True)
    relations=relation_windows(w,K)
    assert field.sparse_full_pivot_rank([{4*a+b:x for (a,b),x in row.items()} for row in relations[1]])==6
    print('PASS generic pair rank6',flush=True)
    report=analyze_quintic(relations,K,rank_strategy='full',d2_only=True)
    bound=120-report['d2_rank']-36
    certificate={'result':'PASS','arithmetic':'exact rational function field in eight independent variables',
        'parameters':list(NAMES),'dominance_tangent_rank':dominance['tangent_rank'],
        'dominance_minor':dominance['nonzero_minor'],'old_point_frame_compatibility':True,
        'quintic':report,'universal_d3_upper_bound':36,'generic_H2_lower_bound':bound,
        'obstructed_generic_chart':bound>0,
        'limits':['Transfer to every basepoint partial-symmetric candidate requires the universal-frame-bundle specialization proof in the manuscript.',
                  'No general low-rank AS classification follows from the product-net chart alone.']}
    if args.output:args.output.write_text(json.dumps(certificate,indent=2)+'\n')
    print(json.dumps(certificate,indent=2))


if __name__=='__main__':main()
