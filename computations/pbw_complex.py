"""Exact linearized quadratic PBW complexes (Python 3 + SymPy).

Generators are 0<1<2<3. Rows: lexicographic relation pairs, then
lexicographic ordered degree-2 monomials; cubic rows similarly use
strictly descending triples and ordered degree-3 monomials.
No floating-point arithmetic or network access is used.
"""
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
import sympy as s

PAIRS=list(combinations(range(4),2))
MON2=list(combinations_with_replacement(range(4),2))
MON3=list(combinations_with_replacement(range(4),3))
TRIPLES=list(combinations(range(4),3))

def add(d,k,v):
    d[k]=d.get(k,0)+v
    if d[k]==0: del d[k]

def matrices(rules):
    # rules[(b,a)] is the ordered quadratic expression for x_b x_a, a<b.
    @lru_cache(None)
    def nf(w):
        for j in range(len(w)-1):
            if w[j]>w[j+1]:
                d={}
                for m,c in rules[w[j:j+2]].items():
                    for n,v in nf(w[:j]+m+w[j+2:]).items(): add(d,n,c*v)
                return d
        return {w:s.Integer(1)}
    @lru_cache(None)
    def lin(w):
        # linear term of normal form; key = (word, relation start 0 or 1, parameter)
        for j in range(len(w)-1):
            if w[j]>w[j+1]:
                d={}
                for m,c in rules[w[j:j+2]].items():
                    for k,v in lin(w[:j]+m+w[j+2:]).items(): add(d,k,c*v)
                r=PAIRS.index(tuple(reversed(w[j:j+2])))
                for l,m in enumerate(MON2):
                    for n,v in nf(w[:j]+m+w[j+2:]).items(): add(d,(n,j,10*r+l),v)
                return d
        return {}
    def forced(w,j):
        d={}
        for m,c in rules[w[j:j+2]].items():
            for k,v in lin(w[:j]+m+w[j+2:]).items(): add(d,k,c*v)
        r=PAIRS.index(tuple(reversed(w[j:j+2])))
        for l,m in enumerate(MON2):
            for n,v in nf(w[:j]+m+w[j+2:]).items(): add(d,(n,j,10*r+l),v)
        return d
    M=[s.zeros(80,60),s.zeros(80,60)]
    for r,t in enumerate(TRIPLES):
        w=tuple(reversed(t));d=forced(w,0)
        for k,v in forced(w,1).items():add(d,k,-v)
        for (n,j,c),v in d.items():M[j][r*20+MON3.index(n),c]=v
    N=[s.zeros(60,16),s.zeros(60,16)]
    for r,(a,b) in enumerate(PAIRS):
        relation={(b,a):s.Integer(1)}
        for m,c in rules[(b,a)].items():add(relation,m,-c)
        for j in range(2):
            for source in range(4):
                for target in range(4):
                    d={}
                    for w,c in relation.items():
                        if w[j]!=source:continue
                        m=w[:j]+(target,)+w[j+1:]
                        for n,v in nf(m).items():add(d,n,c*v)
                    for m,v in d.items():N[j][r*10+MON2.index(m),source*4+target]=v
    assert M[0]*N[0]==s.zeros(80,16)
    assert M[1]*N[1]==s.zeros(80,16)
    assert M[0]*N[1]+M[1]*N[0]==s.zeros(80,16)
    return N,M

def quantum(q):return {(b,a):{(a,b):s.Rational(q[r])} for r,(a,b) in enumerate(PAIRS)}

def clifford0():
    # W=0,Y=1,Z=2,X=3, as in the existing period-2 PBW family at lambda=0.
    return {
      (1,0):{(0,1):-1}, (2,0):{(0,2):-1},
      (2,1):{(1,2):-1,(0,0):1},
      (3,0):{(0,3):-1,(2,2):1},
      (3,1):{(1,3):-1,(2,2):1,(0,0):1},
      (3,2):{(2,3):-1},
    }

def rank(A):
    from sympy.polys.matrices import DomainMatrix
    return DomainMatrix.from_Matrix(A).to_field().rank()



def periodic_matrices(rs,m):
    assert len(rs)==m
    @lru_cache(None)
    def nf(w,start):
        for j in range(len(w)-1):
            pair=w[j:j+2]
            if pair[0]>pair[1]:
                d={}
                for mon,c in rs[(start+j)%m][pair].items():
                    for n,v in nf(w[:j]+mon+w[j+2:],start).items():add(d,n,c*v)
                return d
        return {w:s.Integer(1)}
    @lru_cache(None)
    def lin(w,start):
        for j in range(len(w)-1):
            pair=w[j:j+2]
            if pair[0]>pair[1]:
                d={}
                for mon,c in rs[(start+j)%m][pair].items():
                    for k,v in lin(w[:j]+mon+w[j+2:],start).items():add(d,k,c*v)
                r=PAIRS.index(tuple(reversed(pair)))
                for l,mon in enumerate(MON2):
                    for n,v in nf(w[:j]+mon+w[j+2:],start).items():add(d,(n,60*((start+j)%m)+10*r+l),v)
                return d
        return {}
    def forced(w,start,j):
        d={};pair=w[j:j+2]
        for mon,c in rs[(start+j)%m][pair].items():
            for k,v in lin(w[:j]+mon+w[j+2:],start).items():add(d,k,c*v)
        r=PAIRS.index(tuple(reversed(pair)))
        for l,mon in enumerate(MON2):
            for n,v in nf(w[:j]+mon+w[j+2:],start).items():add(d,(n,60*((start+j)%m)+10*r+l),v)
        return d
    M=s.zeros(80*m,60*m);N=s.zeros(60*m,16*m)
    for start in range(m):
        for r,t in enumerate(TRIPLES):
            w=tuple(reversed(t));d=forced(w,start,0)
            for k,v in forced(w,start,1).items():add(d,k,-v)
            for (n,c),v in d.items():M[start*80+r*20+MON3.index(n),c]=v
        for r,(a,b) in enumerate(PAIRS):
            rel={(b,a):s.Integer(1)}
            for mon,c in rs[start][(b,a)].items():add(rel,mon,-c)
            for j in range(2):
                for source in range(4):
                    for target in range(4):
                        d={}
                        for w,c in rel.items():
                            if w[j]!=source:continue
                            mon=w[:j]+(target,)+w[j+1:]
                            for n,v in nf(mon,start).items():add(d,n,c*v)
                        for mon,v in d.items():N[start*60+r*10+MON2.index(mon),16*((start+j)%m)+source*4+target]=v
    assert M*N==s.zeros(80*m,16*m)
    return N,M

