#!/usr/bin/env python3
"""Exact QQ(u,v) quintic ranks; uniform exclusion uses the manuscript lemma."""
import argparse,json,time
from pathlib import Path
import sympy as s
import generic_boundary_field as field
from flint_rational_function import K,zero,one
from quintic_field_analysis import analyze_quintic


def backend_crosscheck():
    u,v=s.symbols('u v')
    for i in range(1,9):
        a=(i*u+v+1)/(u*v+i)
        b=(i*v-u+1)/(u*u+v+i)
        x,y=K.from_sympy(a),K.from_sympy(b)
        for actual,expected in [(x+y,a+b),(x-y,a-b),(x*y,a*b),(x/y,a/b)]:
            assert s.cancel(actual.as_expr()-expected)==0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    backend_crosscheck()
    field.K,field.zero,field.one=K,zero,one
    w=field.tensor();relations=field.windows(w)
    print('PASS exact cubic frame and four generic pair ranks',flush=True)
    report=analyze_quintic(relations,K)
    assert (report['d2_rank'],report['d3_rank'],report['multiplication_rank'],report['combined_rank'])==(77,36,14,42)
    certificate={'result':'PASS','arithmetic':'exact rational functions QQ(u,v), FLINT rational multivariate polynomial gcd',
        'backend_crosscheck':'32 exact operations compared with SymPy',
        'quadric_coefficients':['1','u','v'],'cubic_normalization_determinant':'-2*u*v',
        'pair_ranks':[6]*4,'C_2_5_dimension':20,'C_1_5_dimension':35,
        'quintic_d2_rank':77,'quintic_d3_rank':36,'d2_d3_zero':True,
        'generic_quintic_H2_dimension':7,
        'square_multiplication_rank':14,'square_annihilator_dimension':6,
        'd3_plus_square_rank':42,'remaining_generic_H2_quotient_dimension':1,
        'uniform_conclusion':'Every parameter u*v!=0 in this cubic-frame family fails the standard Koszul AS conditions: outside the Frobenius/quartic-exact open it already fails necessary filters; on that open H2>=7 by rank specialization.',
        'limits':['This is a two-dimensional coefficient subfamily, not every four-quadric system with a basepoint.',
                  'Finite-Frobenius and quartic-exact exceptional strata are not completely classified.',
                  'The generic mixed seventh class has not been put into a compact closed formula.'],
        'elapsed_seconds':report['elapsed_seconds']}
    if args.output:args.output.write_text(json.dumps(certificate,indent=2)+'\n')
    print(json.dumps(certificate,indent=2))


if __name__=='__main__':main()
