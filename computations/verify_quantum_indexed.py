#!/usr/bin/env python3
"""Full nonzero quantum indexed criterion: finite permutations and certificates."""
from fractions import Fraction as F
from itertools import combinations, permutations

PAIRS=list(combinations(range(4),2))
TRIPLES=list(combinations(range(4),3))


def entry(q,a,b):
    return F(1) if a==b else (F(q[a,b]) if a<b else 1/F(q[b,a]))


def triangle(q,a,b,c):
    return entry(q,a,b)*entry(q,b,c)*entry(q,c,a)


def indexed_witness(q,other):
    assert all(q[p] and other[p] for p in PAIRS)
    for permutation in permutations(range(4)):
        if all(triangle(q,a,b,c)==triangle(other,permutation[a],permutation[b],permutation[c])
               for a,b,c in TRIPLES):
            d=[F(1)]+[entry(q,0,b)/entry(other,permutation[0],permutation[b])
                       for b in (1,2,3)]
            assert all(entry(other,permutation[a],permutation[b])
                       ==entry(q,a,b)*d[a]/d[b] for a,b in PAIRS)
            return permutation,d
    return None


def sample(x,y,z):
    return {(0,1):F(1),(1,2):F(1),(2,3):F(1),
            (0,2):F(x),(1,3):F(y),(0,3):F(z)}


def main():
    fixtures=[sample(1,1,1),sample(1,2,2),sample(1,2,3),sample(2,3,5)]
    assert [sum(triangle(q,*triple)==1 for triple in TRIPLES) for q in fixtures]==[4,2,1,0]
    for q in fixtures:
        for permutation in permutations(range(4)):
            diagonal=[F(2),F(3),F(5),F(7)]
            transformed={}
            for a,b in PAIRS:
                left,right=permutation[a],permutation[b]
                value=entry(q,a,b)*diagonal[a]/diagonal[b]
                transformed[min(left,right),max(left,right)]=value if left<right else 1/value
            witness=indexed_witness(q,transformed)
            assert witness is not None
            sigma,d=witness
            # Verify relation ratios in positive and negative indexed positions.
            for i in (-3,0,1,4):
                for a,b in PAIRS:
                    ratio=entry(q,a,b)*(d[a]**(-i)*d[b]**(-i-1))/(d[b]**(-i)*d[a]**(-i-1))
                    assert ratio==entry(transformed,sigma[a],sigma[b])
        print('PASS support type triangles=',sum(triangle(q,*t)==1 for t in TRIPLES),
              ': all 24 permutations and diagonal indexed certificates')
    assert indexed_witness({pair:F(1) for pair in PAIRS},
                           {pair:F(-1) for pair in PAIRS}) is None
    assert indexed_witness(sample(2,2,1),sample(3,3,1)) is None
    print('PASS rho=6 two-class separation and nonisomorphic equal-rho=8 control')
    import sympy as sym
    p,q,r=sym.symbols('p q r',nonzero=True)
    transformations=[lambda x:(1/x[0],1/x[1],x[0]*x[2]/x[1]),
                     lambda x:(1/x[0],x[2],x[1]),
                     lambda x:(x[1],x[0],1/x[2])]
    original=(p,q,r)
    def equal(left,right):
        return all(sym.cancel(a-b)==0 for a,b in zip(left,right))
    for transform in transformations:
        assert equal(transform(transform(original)),original)
    for i in (0,1):
        a,b=transformations[i:i+2]
        assert equal(a(b(a(original))),b(a(b(original))))
    assert equal(transformations[0](transformations[2](original)),
                 transformations[2](transformations[0](original)))
    for source in fixtures:
        normalized=sample(1,1,1)
        d=[F(1)]+[entry(source,0,b) for b in (1,2,3)]
        normalized={pair:entry(source,*pair)*d[pair[0]]/d[pair[1]] for pair in PAIRS}
        assert all(normalized[0,b]==1 for b in (1,2,3))
        assert [normalized[pair] for pair in ((1,2),(1,3),(2,3))]==[
            triangle(source,*triple) for triple in ((0,1,2),(0,1,3),(0,2,3))]
        assert indexed_witness(source,normalized) is not None
    print('PASS three-torus normalization and symbolic S4 Coxeter relations')
    print('PASS exact rational certificates; geometric necessity proved in manuscript')


if __name__=='__main__':
    main()
