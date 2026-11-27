"""Intrinsic frames, coefficient zeros, and unmarked orbits on the full 112 period-two conic."""
from functools import reduce
from itertools import combinations,permutations,product
from math import prod

from gaussian_curve_field import Gaussian,ONE,convert,h,strip,zero
from verify_parallel_invertible_rectangles import K,WORDS,X,s,smith,template
from verify_rank_fourteen_conic_presentation import conic_data
from verify_third_rectangle_orbits import support_maps


def geometry():
    i = s.I
    S,E = (h-2/h)/2,(h+2/h)/(2*i)
    r = (S*S+1)/2
    l,m = ((r+1)*S+(r-1)*E)/(2*r),((r+1)*S-(r-1)*E)/(2*r)
    A,B,C,D = h*h+2*h-2,h*h-2*h-2,h*h+2*i*h+2,h*h-2*i*h+2
    u,v,z = -A*C/(B*D),-A*D/(B*C),i*(h**4-4)/(4*h*h)
    def same(a,b):
        return s.cancel(a-b,extension=i)==0
    assert same((1+l)/(1-l),u) and same((1+m)/(1-m),v) and same(-S*E,z)
    assert all(same(f.subs(h,i*h),g) for f,g in zip((u,v,z),(v,1/u,-z)))
    assert all(same(f.subs(h,-h),g) for f,g in zip((u,v,z),(1/u,1/v,z)))

    def polynomial(f):
        return s.Poly(f,h,extension=i)
    punctures = h*(h**4+8*h*h+4)
    for sign in (-1,1):
        punctures *= (h*h+sign*2*h-2)*(h*h+sign*2*h+2)*(h**4+sign*4*h**3+8*h*h-sign*8*h+4)
    boundary = polynomial(1)
    for f in (u,v,u+v,u*v+1,z*z-1,z+v,z-u,u*z-1,v*z+1):
        for q in s.fraction(s.cancel(f,extension=i)):
            boundary = boundary.lcm(polynomial(q))
    assert boundary.sqf_part().monic()==polynomial(punctures).monic() and s.degree(punctures,h)==21

    _,W,_,_ = conic_data()
    base = template(112)
    holes = [k for k,f in base.items() if f in X]
    exceptions = {}
    all_factors = set()
    for j,k in enumerate(holes):
        numerator = polynomial(s.fraction(s.cancel(W[k],extension=i))[0])
        f = numerator.sqf_part().monic()
        f = f.exquo(f.gcd(boundary))
        if f.degree():
            assert numerator.gcd(numerator.diff()).gcd(f).degree()==0
            factors = {q.monic().as_expr() for q,_ in f.factor_list()[1]}
            assert len(factors)==2 and all(s.degree(q,h)==3 for q in factors)
            assert not all_factors.intersection(factors)
            all_factors.update(factors)
            exceptions[j] = set(base)-{k}
    assert set(exceptions)=={3,6,8,10} and len(all_factors)==8
    for source in [set(base),*exceptions.values()]:
        weights = s.Matrix([[int(k[p]==j) for p,j in product(K,repeat=2)] for k in sorted(source)])
        diag,_,_ = smith(weights)
        assert [abs(v) for v in diag.diagonal() if v]==[1]*11+[2,2]
        for target in [set(base),*exceptions.values()]:
            if len(source)==len(target):
                for pi in support_maps(source,target):
                    assert {tuple(p[i] for p,i in zip(pi,k)) for k in base}==set(base)


def frames():
    _,raw,units,_ = conic_data()
    coeff = {k:convert(f) for k,f in raw.items()}
    den = ONE
    for x in coeff.values():
        den *= x.d//x.d.gcd(den)
    W = {k:Gaussian(x.a*(den//x.d),x.b*(den//x.d)) for k,x in coeff.items()}
    def sign(seq):
        return (-1)**sum(a>b for a,b in combinations(seq,2))
    G = {}
    for k in WORDS:
        complement = [tuple(j for j in K if j!=i) for i in k]
        value = zero
        for tail in product(*(permutations(c) for c in complement[1:])):
            slots = (complement[0],)+tail
            a,b,c = zip(*slots)
            if a in W and b in W and c in W:
                value += W[a]*W[b]*W[c]*prod(sign((i,)+row) for i,row in zip(k,slots))
        G[k] = value
    for p in K:
        moments = []
        for q in K:
            if p==q:
                continue
            rest = tuple(j for j in K if j not in (p,q))
            def key(g,i,b):
                word = [0]*4
                for j,k in zip((p,q)+rest,(i,i^g,b,b^g)):
                    word[j] = k
                return tuple(word)
            diagonal = [zero]*4
            for g in K:
                block = [[sum((W.get(key(g,i,b),zero)*G[key(g,j,b)] for b in K),zero)
                          for j in K] for i in K]
                for i in K:
                    diagonal[i] += sum((block[i][j]*block[j][i] for j in K),zero)
            moments.append(diagonal)
        for i,j in combinations(K,2):
            gcd = reduce(lambda a,b:a.gcd(b),((row[i]-row[j]).norm() for row in moments))
            assert strip(gcd,units).degree()==0


if __name__=='__main__':
    geometry()
    frames()
    print('PASS: the 22 punctures, 24 simple support defects, all intrinsic frames and their complete normalizations')
