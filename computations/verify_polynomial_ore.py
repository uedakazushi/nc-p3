#!/usr/bin/env python3
"""Symbolic PBW/trace/rank certification for all 18 polynomial Ore coefficients."""
import argparse
import json
from pathlib import Path
from fractions import Fraction
from itertools import combinations,combinations_with_replacement,product
from frobenius_tensor import top_tensor, flatten
from tensor_symmetrization import symmetrize
from functools import lru_cache
import sympy as s
P=list(combinations_with_replacement(range(3),2))
C=s.symbols(' '.join(f'c{i}_{a}{b}' for i in range(3) for a,b in P))
def coefficient(i,a,b):return C[6*i+P.index((min(a,b),max(a,b)))]
def rules(a,b):
    if a<b:return None
    if a==b:
        return [] if a==3 else [((i,3),coefficient(i,a,a)) for i in range(3)]
    ans=[((b,a),s.Integer(-1))]
    if a!=3:ans += [((i,3),coefficient(i,a,b)) for i in range(3)]
    return ans
@lru_cache(None)
def nf(word):
    for pos in range(len(word)-1):
        rs=rules(word[pos],word[pos+1])
        if rs is not None:
            out={}
            for pair,co in rs:
                for normal,value in nf(word[:pos]+pair+word[pos+2:]).items():
                    out[normal]=out.get(normal,0)+co*value
            return {w:s.expand(x) for w,x in out.items() if s.expand(x)!=0}
    return {word:s.Integer(1)}
def ambiguity(word,pos):
    out={}
    for pair,co in rules(word[pos],word[pos+1]):
        for w,value in nf(word[:pos]+pair+word[pos+2:]).items():out[w]=out.get(w,0)+co*value
    return {w:s.expand(x) for w,x in out.items() if s.expand(x)!=0}
def idx(word):
    val=0
    for i in word:val=4*val+i
    return val


def primal_window(coefficients):
    relations=[{(b,a):Fraction(1),(a,b):Fraction(-1)}
               for a in range(3) for b in range(a+1,3)]
    for i in range(3):
        relation={(3,i):Fraction(1),(i,3):Fraction(-1)}
        for k,pair in enumerate(P):
            value=coefficients[6*i+k]
            if value:
                relation[pair]=-Fraction(value)
        relations.append(relation)
    return relations


def primal_rules(a,b):
    if a<=b:
        return None
    answer=[((b,a),s.Integer(1))]
    if a==3:
        answer += [(pair,coefficient(b,*pair)) for pair in P]
    return answer


@lru_cache(None)
def primal_nf(word):
    for pos in range(len(word)-1):
        replacements=primal_rules(word[pos],word[pos+1])
        if replacements is not None:
            result={}
            for pair,value in replacements:
                for normal,co in primal_nf(word[:pos]+pair+word[pos+2:]).items():
                    result[normal]=result.get(normal,0)+value*co
            return {normal:s.expand(co) for normal,co in result.items() if s.expand(co)!=0}
    return {word:s.Integer(1)}


def primal_ambiguity(word,pos):
    result={}
    for pair,value in primal_rules(word[pos],word[pos+1]):
        for normal,co in primal_nf(word[:pos]+pair+word[pos+2:]).items():
            result[normal]=result.get(normal,0)+value*co
    return {normal:s.expand(co) for normal,co in result.items() if s.expand(co)!=0}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    # Verify all dual rewrite relations annihilate the six primal relations
    # as polynomial identities, including the diagonal/off-diagonal factors.
    primal=[]
    for b in range(1,4):
        for a in range(b):
            relation={(b,a):s.Integer(1),(a,b):s.Integer(-1)}
            if b==3:
                for pair in P:
                    relation[pair]=-coefficient(a,*pair)
            primal.append(relation)
    dual=[]
    for b in range(4):
        for a in range(b,4):
            relation={(a,b):s.Integer(1)}
            for pair,value in rules(a,b):
                relation[pair]=relation.get(pair,0)-value
            dual.append(relation)
    assert len(primal)==6 and len(dual)==10
    assert all(s.expand(sum(value*d.get(pair,0) for pair,value in r.items()))==0
               for r in primal for d in dual)
    primal_overlaps=0
    for word in product(range(4),repeat=3):
        if word[0]>word[1]>word[2]:
            assert primal_ambiguity(word,0)==primal_ambiguity(word,1),word
            primal_overlaps+=1
    assert primal_overlaps==4
    dual_overlaps=0
    for word in product(range(4),repeat=3):
        if word[0]>=word[1]>=word[2]:
            assert ambiguity(word,0)==ambiguity(word,1),word
            dual_overlaps+=1
    assert dual_overlaps==20
    tensor={idx(word):nf(word).get((0,1,2,3),s.Integer(0))
            for word in product(range(4),repeat=4)}
    tensor={i:x for i,x in tensor.items() if x!=0}
    assert tensor[27]==1
    assert all(s.Poly(x,*C).total_degree()<=1 for x in tensor.values())
    # The only coefficient-dependent words have no generator 3.
    assert all(not x.free_symbols or 3 not in (i//64,(i//16)%4,(i//4)%4,i%4)
               for i,x in tensor.items())
    def orders(pair):
        return [pair] if pair[0]==pair[1] else [pair,pair[::-1]]
    q=s.Matrix(6,6,lambda i,j:sum(tensor.get(idx((a,b,c,d)),0)
        for a,c in orders(P[i]) for b,d in orders(P[j])))
    assert q.T==-q
    c=lambda i,ab:coefficient(i,int(ab[0]),int(ab[1]))
    expected=[2*c(2,'00'),-2*c(1,'00'),c(2,'01'),-c(1,'01')+c(2,'02'),-c(1,'02'),
        2*c(0,'00')-c(1,'01')-c(2,'02'),2*c(2,'11'),
        c(0,'01')-2*c(1,'11')+c(2,'12'),c(0,'02')-c(1,'12'),
        -c(0,'01')+c(2,'12'),-c(0,'02')-c(1,'12')+2*c(2,'22'),-2*c(1,'22'),
        2*c(0,'11'),c(0,'12'),2*c(0,'22')]
    upper=list(combinations(range(6),2))
    assert [q[pair] for pair in upper]==expected
    linear=s.Matrix(15,18,lambda i,j:q[upper[i]].coeff(C[j]))
    assert linear.rank()==15
    radial=s.zeros(18,3)
    for i in range(3):
        for j in range(3):
            radial[6*i+P.index(tuple(sorted((i,j)))),j]=1
    assert radial.rank()==3 and linear*radial==s.zeros(15,3)
    pivots=linear.rref()[1]
    minor=linear[:,list(pivots)]
    inverse=minor.inv()
    fixtures=[]
    for desired in (0,2,4,6):
        target=s.zeros(6)
        for i in range(0,desired,2):
            target[i,i+1]=1;target[i+1,i]=-1
        values=inverse*s.Matrix([target[pair] for pair in upper])
        coefficients=[s.Integer(0)]*18
        for column,value in zip(pivots,values):
            coefficients[column]=value
        substitution=dict(zip(C,coefficients))
        actual={i:Fraction(x.subs(substitution)) for i,x in tensor.items() if x.subs(substitution)!=0}
        rel=primal_window(coefficients)
        independently=top_tensor([rel]*4)
        scale=independently[27]
        independently={i:x/scale for i,x in independently.items()}
        assert independently==actual
        assert flatten(actual,(0,2))==6+desired
        assert q.subs(substitution).rank()==desired
        alternating,_,_=symmetrize(actual,-1)
        assert (alternating is not None)==(desired==0)
        fixtures.append({'skew_rank':desired,'opposite_rank':6+desired,
                         'coefficients':[str(x) for x in coefficients]})
        print(f'PASS rank stratum: rank Q={desired}, rho={6+desired}',flush=True)
    # A nonzero radial field is the polynomial indexed orbit too.
    coefficients=list(radial*s.Matrix([1,2,-1]))
    substitution=dict(zip(C,coefficients))
    actual={i:Fraction(x.subs(substitution)) for i,x in tensor.items() if x.subs(substitution)!=0}
    assert flatten(actual,(0,2))==6
    assert symmetrize(actual,-1)[0] is not None
    x0,x1,x2,X,n=s.symbols('x0 x1 x2 X n')
    ell=x0+2*x1-x2
    for x in (x0,x1,x2):
        assert s.expand((X-n*ell)*x-x*(X-(n+1)*ell)-ell*x)==0
    certificate={'result':'PASS','arithmetic':'exact rational polynomial identities',
        'primal_cubic_overlaps':primal_overlaps,'dual_cubic_overlaps':dual_overlaps,
        'top_tensor':'linear in all 18 coefficients, verified by dual normal forms',
        'skew_map_rank':15,'kernel_dimension':3,'kernel':'delta_i=x_i*ell',
        'nonzero_minor_columns':list(pivots),'nonzero_minor_determinant':str(minor.det()),
        'rank_strata':fixtures,'nonzero_radial_indexed_witness':'PASS',
        'limits':['No claim of 18-dimensional indexed image or a full moduli component.',
                  'This Ore chart does not cover all standard AS four windows.',
                  'No finite-Frobenius non-Koszul example is produced.']}
    if args.output:
        args.output.write_text(json.dumps(certificate,ensure_ascii=False,indent=2)+'\n')
    print('PASS: all-parameter symbolic PBW/dual identities, skew map/kernel and independent tensor checks')


if __name__=='__main__':
    main()
