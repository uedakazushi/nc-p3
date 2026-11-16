#!/usr/bin/env python3
"""Exact diagonal quantum Ore derivation conditions, gauges and rank controls."""
import argparse
import json
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
import sympy as s
from frobenius_tensor import top_tensor,flatten
from tensor_symmetrization import symmetrize
from verify_quartic_obstruction import relation_matrix,insertion_columns,exact_rank

PAIRS=list(combinations_with_replacement(range(3),2))
CUBIC=list(combinations_with_replacement(range(3),3))


def normal(word,q):
    coefficient=s.Integer(1)
    for i in range(len(word)):
        for j in range(i+1,len(word)):
            if word[i]>word[j]:coefficient*=q[word[j],word[i]]
    return tuple(sorted(word)),coefficient


def derivation_matrix(q):
    result=s.zeros(30,18)
    for row,(a,b) in enumerate(combinations(range(3),2)):
        terms=[(b,False,a,1),(a,True,b,q[b,3]),
               (a,False,b,-q[a,b]),(b,True,a,-q[a,b]*q[a,3])]
        for component,front,letter,scale in terms:
            for k,pair in enumerate(PAIRS):
                word=(letter,)+pair if front else pair+(letter,)
                sorted_word,coefficient=normal(word,q)
                result[10*row+CUBIC.index(sorted_word),6*component+k]+=scale*coefficient
    return result.applyfunc(s.expand)


def relations(q,c):
    result=[{(b,a):F(1),(a,b):-F(q[a,b])} for a,b in combinations(range(3),2)]
    for a in range(3):
        r={(3,a):F(1),(a,3):-F(q[a,3])}
        for k,pair in enumerate(PAIRS):
            if c[6*a+k]:r[pair]=-F(c[6*a+k])
        result.append(r)
    return result


def inner_coefficients(q,ell):
    result=s.zeros(18,1)
    for i in range(3):
        for j in range(3):
            pair,scale=normal((j,i),q)
            result[6*i+PAIRS.index(pair)]+=ell[j]*scale
            pair,scale=normal((i,j),q)
            result[6*i+PAIRS.index(pair)]-=q[i,3]*ell[j]*scale
    return result


def actual_tensor(q,c):
    rel=relations(q,c)
    w=top_tensor([rel]*4)
    return w


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    pairs=list(combinations(range(4),2))
    variables=s.symbols('q01 q02 q03 q12 q13 q23',nonzero=True)
    symbolic=dict(zip(pairs,variables))
    matrix=derivation_matrix(symbolic)
    # Universal inner sigma-derivations provide independent identity controls.
    inner=s.Matrix.hstack(*(inner_coefficients(symbolic,[int(j==i) for j in range(3)])
                            for i in range(3)))
    assert all(s.expand(x)==0 for x in matrix*inner)
    fixtures=[]
    generic={pair:s.Rational(2+pair[0]+2*pair[1]) for pair in pairs}
    cases=[('generic_inner',generic,inner_coefficients(generic,[1,2,3]),12,False,False)]
    for name,sign,diagonal,monomial,rho,issym,isalt in [
        ('commuting_resonance',1,[1,2,3,6],(1,2),8,False,False),
        ('anticommuting_central',-1,[1,2,3,4],(1,1),6,True,False),
        ('anticommuting_outer',-1,[1,2,3,-6],(1,2),8,False,False)]:
        q={(a,b):s.Rational(sign*diagonal[b],diagonal[a]) for a,b in pairs}
        coefficients=s.zeros(18,1);coefficients[PAIRS.index(monomial)]=1
        cases.append((name,q,coefficients,rho,issym,isalt))
    for name,q,c,expected,symmetric,alternating in cases:
        condition=derivation_matrix(q)
        assert condition*c==s.zeros(30,1),(name,condition*c)
        w=actual_tensor(q,c)
        background=actual_tensor(q,[0]*18)
        rho=flatten(w,(0,2));floor=flatten(background,(0,2))
        assert rho==expected and floor<=rho
        assert (symmetrize(w)[0] is not None)==symmetric
        assert (symmetrize(w,-1)[0] is not None)==alternating
        if name!='generic_inner':
            sign=1 if name=='commuting_resonance' else -1
            diagonal=[s.Integer(1)]+[q[0,a]/sign for a in (1,2,3)]
            assert all(s.simplify(q[a,b]*diagonal[a]/diagonal[b]-sign)==0 for a,b in pairs)
            uniform={pair:s.Integer(sign) for pair in pairs}
            for n in (-2,-1,0,1,2,3):
                f=s.zeros(18,1)
                for a in range(3):
                    for k,(b,d) in enumerate(PAIRS):
                        f[6*a+k]=c[6*a+k]*(diagonal[a]/diagonal[d])*(diagonal[3]*diagonal[a]/(diagonal[b]*diagonal[d]))**n
                maps=[s.diag(*(value**(-j) for value in diagonal)) for j in (n,n+1)]
                target=relation_matrix(relations(uniform,f))
                annihilator=s.Matrix.hstack(*target.T.nullspace()).T
                image=s.kronecker_product(*maps)*relation_matrix(relations(q,c))
                assert annihilator*image==s.zeros(10,6)
        fixtures.append({'name':name,'derivation_condition_rank':condition.rank(),
                         'background_rho':floor,'rho':rho,'symmetric':symmetric,
                         'alternating':alternating})
    bad=s.zeros(18,1);bad[0]=1
    assert derivation_matrix(generic)*bad!=s.zeros(30,1)
    bad_dim=64-exact_rank(insertion_columns([relations(generic,bad)]*4,0,3))
    assert bad_dim<20
    certificate={'result':'PASS','arithmetic':'exact rational/Laurent polynomial identities',
        'universal_derivation_matrix_shape':[30,18],
        'universal_inner_sigma_derivation_identity':'PASS',
        'rank_and_gauge_fixtures':fixtures,'nonregular_control_cubic_dimension':bad_dim,
        'limits':['Global AS necessity/sufficiency and specialization inequality use the manuscript proof.',
                  'No claim that the quantum Ore chart covers all standard AS objects.']}
    if args.output:args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
