#!/usr/bin/env python3
"""Symbolic quantum top tensor, all opposite-rank strata, diagonal normal forms."""
from itertools import combinations, permutations
import sympy as s


def main():
    pairs=list(combinations(range(4),2))
    q={pair:s.Symbol('q'+''.join(map(str,pair)),nonzero=True) for pair in pairs}
    f=[q[0,2]*q[1,3]-q[0,3]*q[1,2],
       q[0,1]*q[1,2]*q[2,3]-q[0,3],
       q[0,1]*q[1,3]-q[0,2]*q[2,3]]
    def coefficient(word):
        positions={x:i for i,x in enumerate(word)}
        return s.prod(-1 if positions[a]>positions[b] else val
                      for (a,b),val in q.items())
    # Direct membership in every adjacent relation for all 24 coefficients.
    for word in permutations(range(4)):
        for i in range(3):
            a,b=word[i:i+2]
            if a<b:
                swapped=word[:i]+(b,a)+word[i+2:]
                assert s.expand(coefficient(word)+q[a,b]*coefficient(swapped))==0
    blocks=[]
    for k,pair in enumerate(pairs):
        rest=tuple(x for x in range(4) if x not in pair)
        rows=[pair,pair[::-1]];cols=[rest,rest[::-1]]
        block=s.Matrix(2,2,lambda i,j:coefficient(
            (rows[i][0],cols[j][0],rows[i][1],cols[j][1])))
        det=s.factor(block.det());factor=f[[0,1,2,2,1,0][k]]
        unit=s.cancel(det/factor)
        num,den=s.fraction(unit)
        assert den==1 and len(s.Poly(num,*q.values()).terms())==1
        blocks.append(block)
        print(f'PASS block {pair}: det={det}')
    for triple,expected in [((1,1,1),6),((2,2,1),8),((2,3,1),10),((2,3,5),12)]:
        values={q[0,1]:1,q[1,2]:1,q[2,3]:1,
                q[0,2]:triple[0],q[1,3]:triple[1],q[0,3]:triple[2]}
        rank=sum(b.subs(values).rank() for b in blocks)
        assert rank==expected
        assert rank==6+2*sum(x.subs(values)!=0 for x in f)
        print(f'PASS rank stratum {expected}: normalized sample={triple}')
    a,b,c=s.symbols('a b c',nonzero=True)
    for sign in (1,-1):
        values={(0,1):a,(1,2):b,(2,3):c,(0,2):sign*a*b,
                (1,3):sign*b*c,(0,3):a*b*c}
        diagonal=[1,a/sign,a*b,a*b*c/sign]
        assert all(s.simplify(values[i,j]*diagonal[i]/diagonal[j]-sign)==0
                   for i,j in pairs)
        assert all(s.expand(x.subs({q[k]:v for k,v in values.items()}))==0 for x in f)
        print(f'PASS rho=6 diagonal indexed normal form s={sign}')
    print('PASS: exact symbolic identities; no classification outside the quantum subfamily')


if __name__=='__main__':
    main()
