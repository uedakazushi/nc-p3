#!/usr/bin/env python3
"""Exact one-return bimodule Tor regression and optional generic rank."""
import argparse,json,time
from pathlib import Path
import sympy as s
import generic_boundary_field as field
from flint_rational_function import rational_function_field
from product_net_boundary import NAMES,tensor,terms,relation_windows
from trivial_extension_bitor import one_return_complex
from quintic_field_analysis import combine


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field',choices=('point','ci','two','eight'),default='point')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--composition-check',choices=('direct','structural'),default='direct')
    args=parser.parse_args();start=time.monotonic()
    K=rational_function_field(NAMES if args.field=='eight' else ('u','v') if args.field=='two' else ('z',))
    class Coefficients:
        zero,one=K.zero,K.one
        @staticmethod
        def from_sympy(e):
            subs=dict(zip(s.symbols(NAMES),(0,0,0,1,1,0,0,0)))
            if args.field=='eight':subs={}
            if args.field=='two':subs={k:v for k,v in subs.items() if str(k) not in ('u','v')}
            return K.from_sympy(s.sympify(e).subs(subs))
    coefficients=Coefficients();field.K,field.zero,field.one=K,K.zero,K.one
    if args.field=='ci':
        from itertools import combinations_with_replacement
        from verify_nonkoszul_boundary import tensor_at
        w={i:K.from_sympy(v) for i,v in tensor_at(1).items()}
        quadrics=[{(i,i):K.one,(i,3):-K.one} for i in range(3)]
        quadrics.append({(i,3):K.one for i in range(4)})
        pairs=list(combinations_with_replacement(range(4),2))
        matrix=s.Matrix([[q.get(pair,K.zero).as_expr() for pair in pairs] for q in quadrics])
        r0=[]
        for vector in matrix.nullspace():
            row={}
            for pair,value in zip(pairs,vector):
                if value:row[pair]=row[pair[::-1]]=K.from_sympy(value)
            r0.append(row)
    else:
        w=tensor(coefficients);r0=relation_windows(w,coefficients)[0];quadrics=terms(coefficients)
    differentials,bases=one_return_complex(w,r0,quadrics,K)
    print('PASS one-return complex dimensions8,144,480,400',flush=True)
    if args.composition_check=='direct':
        assert all(not combine(differentials[2],v,K.zero) for v in differentials[3])
        assert all(not combine(differentials[3],v,K.zero) for v in differentials[4])
        print('PASS directly checked zero compositions',flush=True)
    else:print('Using proved two-sided Koszul differential; direct composition check omitted',flush=True)
    ranks={h:field.sparse_full_pivot_rank(differentials[h],lambda r:print('rank progress',h,r,flush=True)) for h in (2,3,4)}
    homology=[8-ranks[2],144-ranks[2]-ranks[3],480-ranks[3]-ranks[4],400-ranks[4]]
    if args.field in ('point','two'):assert homology==[0,0,7,63]
    if args.field=='ci':assert homology==[0,0,0,56]
    report={'result':'PASS','field':args.field,'complex_dimensions':[8,144,480,400],
            'differential_ranks':ranks,'homology_dimensions':homology,
            'composition_check':args.composition_check,'direct_zero_compositions_checked':args.composition_check=='direct','nonlinear_one_return_Tor_dimension':sum(homology[:3]),
            'limits':['Only the indicated internal-degree-five vertex summand is calculated.',
                      'No eight-variable result is inferred from a rational point or two-variable slice.'],
            'elapsed_seconds':round(time.monotonic()-start,3)}
    if args.output:args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
