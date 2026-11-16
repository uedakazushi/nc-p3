#!/usr/bin/env python3
"""All-variable PBW, normalization linear algebra and rank for signed Ore windows."""
import argparse
import json
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from itertools import product,combinations,combinations_with_replacement
import sympy as s
from frobenius_tensor import top_tensor,flatten
from tensor_symmetrization import symmetrize
from verify_quartic_obstruction import insertion_columns,exact_rank
P=list(combinations_with_replacement(range(3),2))
H=s.symbols('h0:9');G=s.symbols('g0:3')


def primal_rules(a,b,start):
    if a<=b:return None
    if a<3:return [((b,a),s.Integer(-1))]
    result=[((b,3),s.Integer(-1))]
    result += [((j,j),H[3*b+j]) for j in range(3)]
    result.append((tuple(j for j in range(3) if j!=b),(-1)**start*G[b]))
    return result


def dual_rules(a,b,start):
    if a<b:return None
    if a==b:
        return [] if a==3 else [((j,3),-H[3*j+a]) for j in range(3)]
    if a==3:return [((b,3),s.Integer(1))]
    other=3-a-b
    return [((b,a),s.Integer(1)),((other,3),(-1)**start*G[other])]


def normalizer(rule):
    @lru_cache(None)
    def nf(start,word):
        for pos in range(len(word)-1):
            replacements=rule(word[pos],word[pos+1],(start+pos)%2)
            if replacements is not None:
                result={}
                for pair,value in replacements:
                    for normal,co in nf(start,word[:pos]+pair+word[pos+2:]).items():
                        result[normal]=result.get(normal,0)+value*co
                return {normal:s.expand(co) for normal,co in result.items() if s.expand(co)!=0}
        return {word:s.Integer(1)}
    return nf


def ambiguity(rule,nf,start,word,pos):
    result={}
    for pair,value in rule(word[pos],word[pos+1],(start+pos)%2):
        for normal,co in nf(start,word[:pos]+pair+word[pos+2:]).items():
            result[normal]=result.get(normal,0)+value*co
    return {normal:s.expand(co) for normal,co in result.items() if s.expand(co)!=0}


def anti_nf(word):
    sign=(-1)**sum(word[i]>word[j] for i in range(len(word)) for j in range(i+1,len(word)))
    return tuple(sorted(word)),sign


def adjacent_matrix():
    triples=list(combinations_with_replacement(range(3),3))
    result=s.zeros(30,36)
    for row,(a,b) in enumerate(combinations(range(3),2)):
        for side,component,front,last,sign in [(0,b,False,a,1),(0,a,False,b,1),
                                              (1,a,True,b,-1),(1,b,True,a,-1)]:
            for k,pair in enumerate(P):
                normal,co=anti_nf((last,)+pair if front else pair+(last,))
                result[10*row+triples.index(normal),18*side+6*component+k]+=sign*co
    return result


def window(h,g,start):
    relations=[{(b,a):F(1),(a,b):F(1)} for a in range(3) for b in range(a+1,3)]
    for i in range(3):
        r={(3,i):F(1),(i,3):F(1)}
        for j in range(3):
            if h[3*i+j]:r[j,j]=-F(h[3*i+j])
        if g[i]:r[tuple(j for j in range(3) if j!=i)]=-((-1)**start)*F(g[i])
        relations.append(r)
    return relations


def index(word):
    result=0
    for a in word:result=4*result+a
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    primal_nf=normalizer(primal_rules);dual_nf=normalizer(dual_rules)
    counts=[]
    for rule,nf,weak in [(primal_rules,primal_nf,False),(dual_rules,dual_nf,True)]:
        count=0
        for start in (0,1):
            for word in product(range(4),repeat=3):
                if ((word[0]>=word[1]>=word[2]) if weak else (word[0]>word[1]>word[2])):
                    assert ambiguity(rule,nf,start,word,0)==ambiguity(rule,nf,start,word,1),(start,word)
                    count+=1
        counts.append(count)
    assert counts==[8,40]
    # Verify degree-two orthogonality at both starting parities.
    for start in (0,1):
        primal=[];dual=[]
        for a in range(4):
            for b in range(a):
                r={(a,b):s.Integer(1)}
                for pair,value in primal_rules(a,b,start):r[pair]=r.get(pair,0)-value
                primal.append(r)
            for b in range(a+1):
                r={(a,b):s.Integer(1)}
                for pair,value in dual_rules(a,b,start):r[pair]=r.get(pair,0)-value
                dual.append(r)
        assert all(s.expand(sum(value*d.get(pair,0) for pair,value in r.items()))==0 for r in primal for d in dual)
    tensor={index(word):dual_nf(0,word).get((0,1,2,3),s.Integer(0)) for word in product(range(4),repeat=4)}
    tensor={i:x for i,x in tensor.items() if x!=0}
    assert tensor[27]==1
    assert all(s.Poly(value,*(H+G)).total_degree()<=1 for value in tensor.values())
    # Quotient by off-diagonal symmetric tensors: diagonals and differences.
    domain=[{(i,i):1} for i in range(3)]+[{(a,b):1,(b,a):-1} for a,b in combinations(range(3),2)]
    q=s.Matrix(6,6,lambda i,j:sum(ci*cj*tensor.get(index((a,b,c,d)),0)
        for (a,c),ci in domain[i].items() for (b,d),cj in domain[j].items()))
    expected=s.zeros(6)
    expected[0,1]=-G[2];expected[0,2]=-G[1];expected[1,2]=-G[0]
    expected=expected-expected.T
    assert q==expected,(q,expected)
    adjacent=adjacent_matrix()
    assert adjacent.rank()==18
    pair_parameter=s.zeros(36,21)
    for n in range(2):
        for i in range(3):
            for j in range(3):pair_parameter[18*n+6*i+P.index((j,j)),3*i+j]=1
            pair=tuple(j for j in range(3) if j!=i)
            pair_parameter[18*n+6*i+P.index(pair),9+i]=(-1)**n
            for j in range(3):
                pair,sign=anti_nf((j,i))
                pair_parameter[18*n+6*i+P.index(pair),12+3*n+j]+=sign
                pair,sign=anti_nf((i,j))
                pair_parameter[18*n+6*i+P.index(pair),12+3*(n+1)+j]+=sign
    assert adjacent*pair_parameter==s.zeros(30,21) and pair_parameter.rank()==18
    cycle=s.zeros(120,72)
    for n in range(4):
        for i in range(30):
            for j in range(18):
                cycle[30*n+i,18*n+j]+=adjacent[i,j]
                cycle[30*n+i,18*((n+1)%4)+j]+=adjacent[i,18+j]
    assert cycle.rank()==51
    # Full four-window regular parameterization: H, outer G and X shifts.
    parameter=s.zeros(72,24)
    for n in range(4):
        for i in range(3):
            for j in range(3):parameter[18*n+6*i+P.index((j,j)),3*i+j]=1
            pair=tuple(j for j in range(3) if j!=i)
            parameter[18*n+6*i+P.index(pair),9+i]=(-1)**n
            for j in range(3):
                pair,sign=anti_nf((j,i))
                parameter[18*n+6*i+P.index(pair),12+3*n+j]+=sign
                pair,sign=anti_nf((i,j))
                parameter[18*n+6*i+P.index(pair),12+3*((n+1)%4)+j]+=sign
    assert cycle*parameter==s.zeros(120,24) and parameter.rank()==21
    for gv in ((0,0,0),(1,0,0),(1,2,3)):
        hv=tuple(range(1,10));substitution=dict(zip(H+G,hv+gv))
        actual={i:F(x.subs(substitution)) for i,x in tensor.items() if x.subs(substitution)!=0}
        windows=[window(hv,gv,n) for n in range(4)]
        independent=top_tensor(windows);scale=independent[27]
        assert {i:x/scale for i,x in independent.items()}==actual
        assert flatten(actual,(0,2))==(6 if not any(gv) else 8)
        assert (symmetrize(actual)[0] is not None)==(not any(gv))
        assert all(64-exact_rank(insertion_columns(windows,n,3))==20 for n in range(4))
    certificate={'result':'PASS','arithmetic':'exact rational polynomial/linear algebra',
        'primal_overlaps':counts[0],'dual_overlaps':counts[1],
        'adjacent_rank':18,'complete_adjacent_parameterization_rank':18,'four_window_condition_rank':51,
        'regular_four_window_coefficient_dimension':21,
        'normal_form_parameter_map_rank':21,'normal_form_parameter_count':24,
        'gauge_kernel_dimension':3,'top_tensor':'linear in H and G',
        'opposite_rank':'6 if outer G=0, otherwise 8',
        'low_rank':'exactly the symmetric Clifford subchart',
        'limits':['21 is a presentation coefficient dimension, not indexed moduli dimension.',
                  'No assertion that every nonzero outer G has minimal period two.']}
    if args.output:args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(certificate,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
