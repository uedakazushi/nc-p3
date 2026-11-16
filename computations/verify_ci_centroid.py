#!/usr/bin/env python3
"""Jordan coverage for scalar CI centroids and anchored frame recovery."""
import argparse
import json
from pathlib import Path
from itertools import product,combinations,combinations_with_replacement
from fractions import Fraction as F
import sympy as s
from frobenius_tensor import top_tensor,normalized
from tensor_symmetrization import symmetrizer_space,act_factor
from pbw_complex import PAIRS,clifford0


def partitions(n,maximum=None):
    if n==0:
        yield ();return
    for first in range(min(n,n if maximum is None else maximum),0,-1):
        for tail in partitions(n-first,first):yield (first,)+tail


def jordan(blocks,eigenvalue):
    return s.diag(*(s.Matrix(size,size,lambda i,j:eigenvalue if i==j else int(j==i+1)) for size in blocks))


def minor_forms(matrix,variables):
    images=list(matrix.T*s.Matrix(variables))
    return [s.expand(variables[a]*images[b]-variables[b]*images[a]) for a,b in combinations(range(4),2)]


def from_rules(rules):
    rel=[]
    for a,b in PAIRS:
        r={(b,a):F(1)}
        for pair,value in rules[b,a].items():r[pair]=-F(value)
        rel.append(r)
    return top_tensor([rel]*4)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    x=s.symbols('x0:4');monomials=list(combinations_with_replacement(range(4),2))
    profiles=[];low=[]
    for multiplicities in partitions(4):
        for blocks in product(*(list(partitions(size)) for size in multiplicities)):
            matrix=s.diag(*(jordan(block,i) for i,block in enumerate(blocks)))
            if matrix==s.zeros(4):continue
            forms=minor_forms(matrix,x)
            coefficients=s.Matrix(10,6,lambda i,j:s.Poly(forms[j],*x).coeff_monomial(x[monomials[i][0]]*x[monomials[i][1]]))
            rank=coefficients.rank()
            profile={'multiplicities':list(multiplicities),'Jordan_blocks':[list(b) for b in blocks],'minor_span_rank':rank}
            profiles.append(profile)
            if rank<=4:
                eigen_spaces=[]
                for eigenvalue in range(len(blocks)):
                    kernel=(matrix.T-eigenvalue*s.eye(4)).nullspace()
                    if kernel:eigen_spaces.append(s.Matrix.hstack(*kernel))
                subspace=max(eigen_spaces,key=lambda m:m.cols)
                minimum=3 if rank==3 else 2
                assert subspace.cols>=minimum
                y=s.symbols(f'y0:{subspace.cols}')
                substitution=dict(zip(x,list(subspace*s.Matrix(y))))
                assert all(s.expand(f.subs(substitution,simultaneous=True))==0 for f in forms)
                profile['common_eigenspace_vector_dimension']=subspace.cols
                low.append(profile)
    assert len(profiles)==14
    assert sorted((p['minor_span_rank'],p['common_eigenspace_vector_dimension']) for p in low)==[(3,3),(3,3),(4,2),(4,2)]
    pure={(b,a):{(a,b):-1} for a,b in PAIRS}
    perturbed={key:dict(value) for key,value in pure.items()};perturbed[1,0][3,3]=1
    frames=[]
    for j in range(3):
        frame=s.eye(4);frame[j,j+1]=j+2;frame[(j+2)%4,(j+2)%4]=j+2
        frames.append(frame)
    for name,rules in [('pure',pure),('regular_specialization_generic',perturbed),('C0',clifford0())]:
        w=from_rules(rules)
        assert all(len(symmetrizer_space(w,j))==1 for j in (1,2,3))
        changed=w
        for j,frame in enumerate(frames,1):changed=act_factor(changed,frame,j)
        restored=changed
        for j,frame in enumerate(frames,1):
            space=symmetrizer_space(changed,j)
            assert len(space)==1
            witness=space[0]
            scalar=(witness*frame)[0,0]
            assert scalar and witness*frame==scalar*s.eye(4)
            restored=act_factor(restored,witness,j)
        assert normalized(restored)==normalized(w)
    certificate={'result':'PASS','arithmetic':'exact rational polynomial/linear algebra',
        'nonscalar_Jordan_profiles':len(profiles),'low_minor_rank_profiles':low,
        'anchored_frame_recovery':['pure','regular_specialization_generic','C0'],
        'presentation_product':'PGL4^3 x U, dimension 69',
        'limits':['The Jordan/basepoint-free proof, not the fixtures alone, establishes scalar centroid for every U point.',
                  'Smoothness of the symmetric image does not prove closedness in the finite-Frobenius locus or ambient smoothness at every point.']}
    if args.output:args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
