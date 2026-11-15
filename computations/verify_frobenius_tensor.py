#!/usr/bin/env python3
"""Exact top-tensor reconstruction/rank regression; not an all-AS-locus test."""
import argparse
import json
from itertools import permutations
from fractions import Fraction as F
from pathlib import Path
from frobenius_tensor import (digits, exact_rank, flatten, insertion_columns,
    nullspace_sparse, orthogonal_windows, reconstructed_windows, rotate,
    top_tensor)
from verify_quartic_obstruction import seven_window
from tensor_symmetrization import symmetrize


def nine_window(sign):
    # u=1,v=2,r=3,s=4,a=1,b=2,c=3,t=2,h=3; 1+th=7.
    return [
        {(1, 0): 1, (0, 1): 1, (2, 2): -2},
        {(2, 0): 1, (0, 2): 1, (1, 1): -3},
        {(2, 1): 1, (1, 2): 1, (0, 0): -1},
        {(3, 0): 1, (0, 3): 1, (1, 1): -1-9*sign,
         (2, 2): -2-6*sign, (1, 2): -sign, (0, 0): F(sign, 2),
         (0, 1): 6*sign, (0, 2): 6*sign},
        {(3, 1): 1, (1, 3): 1, (0, 0): -1, (2, 2): -1+sign,
         (1, 1): -3+3*sign, (0, 1): -sign, (0, 2): -2*sign},
        {(3, 2): 1, (2, 3): 1, (1, 1): -4,
         (0, 1): -3*sign, (2, 2): 3*sign},
    ]


def span_rank(relations):
    return exact_rank(({4*a+b: F(x) for (a,b),x in r.items() if x}
                       for r in relations))


def shear(tensor):
    result = dict(tensor)
    # Independent invertible elementary maps on each tensor factor.
    for factor in range(4):
        transformed = {}
        for index, value in result.items():
            word = list(digits(index))
            transformed[index] = transformed.get(index, F(0)) + value
            if word[factor] == factor:
                word[factor] = (factor + 1) % 4
                target = 0
                for letter in word:
                    target = 4 * target + letter
                transformed[target] = transformed.get(target, F(0)) + 2*value
        result = {i: x for i,x in transformed.items() if x}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    poly=[{(b,a):1,(a,b):-1} for a in range(4) for b in range(a+1,4)]
    quantum=[{(b,a):1,(a,b):-F(2+a+2*b)}
             for a in range(4) for b in range(a+1,4)]
    anti=[{(b,a):1,(a,b):1} for a in range(4) for b in range(a+1,4)]
    symmetric_ci=[dict(r) for r in anti]
    symmetric_ci[0][(3,3)]=-1
    samples={'polynomial':[poly]*4,'anticommutative':[anti]*4,
             'symmetric_ci':[symmetric_ci]*4,'quantum':[quantum]*4,
             'seven_signed':[seven_window((-1)**i) for i in range(4)],
             'nine_signed':[nine_window((-1)**i) for i in range(4)]}
    results={}
    for name, windows in samples.items():
        w=top_tensor(windows)
        singles=[flatten(w,(i,)) for i in range(4)]
        adjacent=[flatten(w,(i,(i+1)%4)) for i in range(4)]
        rho=flatten(w,(0,2))
        assert singles==[4]*4 and adjacent==[6]*4
        if name in ('anticommutative','symmetric_ci'):
            for index, value in w.items():
                for word in set(permutations(digits(index))):
                    permuted = 0
                    for letter in word:
                        permuted = 4*permuted+letter
                    assert w.get(permuted, F(0)) == value
        recovered=reconstructed_windows(w)
        # These are the original first three windows. The fourth is the
        # canonical completion, which need not equal a noncanonical input R3.
        for i in range(3):
            assert span_rank(windows[i]+recovered[i])==6
        dual=orthogonal_windows(recovered)
        dims=[]
        for start in range(4):
            local=[len(nullspace_sparse(insertion_columns(dual,start,n),4**n))
                   for n in (3,4)]
            assert local==[4,1], (name,start,local)
            assert top_tensor(recovered,start)==rotate(w,start)
            dims.append(local)
        changed=shear(w)
        assert [flatten(changed,(i,)) for i in range(4)]==singles
        assert [flatten(changed,(i,(i+1)%4)) for i in range(4)]==adjacent
        assert flatten(changed,(0,2))==rho
        symmetric, spaces, witnesses = symmetrize(w)
        expected = name in ('anticommutative', 'symmetric_ci')
        assert (symmetric is not None) == expected
        changed_symmetric, _, changed_witnesses = symmetrize(changed)
        assert (changed_symmetric is not None) == expected
        membership = {'symmetric_chart': expected,
                      'linear_space_dimensions': [len(space) for space in spaces],
                      'basis_shear_invariance': 'PASS'}
        if expected:
            membership['witness_determinants'] = [str(matrix.det()) for matrix in witnesses]
        results[name]={'single_ranks':singles,'adjacent_pair_ranks':adjacent,
                       'opposite_rank':rho,'intersection_dimensions_3_4':dims,
                       'canonical_rotation':'PASS','basis_shear_invariance':'PASS',
                       'linear_symmetrization':membership}
        print(f'PASS {name}: rho={rho}, reconstruction and canonical rotations',flush=True)
    assert [results[n]['opposite_rank'] for n in results]==[6,6,6,12,8,8]
    # Opposite rank one with all singleton ranks four has adjacent rank sixteen.
    rank_one={64*a+16*b+4*a+b:F(1) for a in range(4) for b in range(4)}
    assert flatten(rank_one,(0,2))==1
    assert [flatten(rank_one,(i,)) for i in range(4)]==[4]*4
    assert flatten(rank_one,(0,1))==16
    certificate={'result':'PASS','arithmetic':'exact rational',
        'samples':results,'rank_one_exclusion_control':'PASS',
        'scope':'Finite tensor/Frobenius reconstruction and sampled ranks; no general Koszul certification or new family completeness/novelty claim.'}
    if args.output:
        args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
