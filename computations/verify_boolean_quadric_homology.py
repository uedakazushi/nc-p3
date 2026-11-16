#!/usr/bin/env python3
"""Exact generic quintic homology for the full one-point Boolean quadric open."""
import argparse,json
from pathlib import Path
import generic_boundary_field as field
from flint_rational_function import rational_function_field
from boolean_quadric_frame import polynomial_tensor,relation_windows,data
from quintic_field_analysis import analyze_quintic


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    K=rational_function_field(('u','v','a','b','c'))
    field.K,field.zero,field.one=K,K.zero,K.one
    w=polynomial_tensor(K)
    _,_,quadric=data(K)
    for j,x in ((j,x) for j in range(4) for x in range(4)):
        def coefficient(a,b,c):return w.get(64*a+16*b+4*j+c,K.zero)
        assert sum((q*coefficient(a,b,x) for (a,b),q in quadric.items()),K.zero)==K.zero
    relations=relation_windows(w,K,clear_normalization=False)
    assert len(field.echelon([{4*a+b:x for (a,b),x in row.items()} for row in relations[1]]))==6
    print('PASS exact polynomial cubic frame and generic relation basis',flush=True)
    report=analyze_quintic(relations,K)
    certificate={'result':'PASS','arithmetic':'exact QQ(u,v,a,b,c), no sampling or interpolation',
        'point_values':['1','u','v','a','b','c','a+b+c-1-u-v'],
        'domain':'all seven listed values nonzero',
        'frame':'global evaluation frame rescaled by face products, which are units on the domain',
        'generic_frobenius_quartic_exact':'nonempty open established by the checked old point (1,1,2,2,2)',
        'quintic':report,
        'uniform_AS_exclusion':report['H2_dimension']>0,
        'specialization_argument':'On the Frobenius/quartic-exact open generic ranks bound every specialization; outside it necessary AS conditions already fail.',
        'limits':['The first three quadrics remain the fixed Boolean complete intersection.',
                  'All special finite-Frobenius strata and full indexed isomorphisms are not classified.']}
    if args.output:args.output.write_text(json.dumps(certificate,indent=2)+'\n')
    print(json.dumps(certificate,indent=2))


if __name__=='__main__':main()
