"""Distinguish both rank-twelve curves and all individual unmarked fibers."""
from itertools import combinations, permutations, product

import sympy as s

from verify_transverse_rank_sixteen_isomorphism import K, adjacent, opposite, universal


def alternating_gradient(w):
    def sign(t):
        return (-1)**sum(a>b for a,b in combinations(t,2))
    g = {}
    for key in w:
        lists = [tuple(j for j in K if j!=i) for i in key]
        total = 0
        for rest in product(*(permutations(lst) for lst in lists[1:])):
            slots = [lists[0],*rest]
            coefficient = s.prod(w.get(t,0) for t in zip(*slots))
            if coefficient:
                total += coefficient*s.prod(sign((i,)+seq) for i,seq in zip(key,slots))
        g[key]=s.factor(total)
    return g


def radical_determinant(w):
    z = s.symbols('z:4')
    ker = s.zeros(4)
    n = 0
    for b in K:
        for col in opposite(w,b).T.nullspace():
            for i in K:
                ker[i,i^b]+=z[n]*col[i]
            n += 1
    assert n==4
    return s.factor(ker.det())


def first_curve(q,z,w,g):
    assert s.factor(radical_determinant(w)+z[2]**2*z[3]**2)==0
    X = q*q
    cap = 16*(X+1)
    L = [adjacent(w,b)*adjacent(g,b).T/cap for b in K]
    fs = [(16*X**2+40*X+11)/(8*(X+1)**2),
          (64*X**3+336*X**2+252*X+31)/(32*(X+1)**3)]
    for i,f in enumerate(fs,2):
        assert s.factor(sum((m**i).trace() for m in L)/2-f)==0
    values = [s.factor(sum((m**2)[j,j] for m in L)) for j in K]
    assert values[0]==values[2] and values[1]==values[3]
    third = (4*X-1)**2/(256*(X+1)**4)
    assert s.factor((values[0]-values[1])**2-third)==0
    a,b,inv = s.symbols('a b inv')
    polys = []
    for f in [*fs,third]:
        ff = s.factor(f.subs(q,s.sqrt(a)))
        polys.append(s.fraction(s.factor((ff-ff.subs(a,b))/(a-b)))[0])
    assert list(s.groebner(polys+[inv*(a-b)*(a+1)*(b+1)-1],inv,a,b))==[1]


def second_curve(p,z,w,g):
    X = p**4
    rad = (z[1]**2-z[2]*z[3])*(X*(z[2]+z[3])**2+4*p**2*z[0]**2-(z[2]-z[3])**2)/(4*p**2)
    assert s.factor(radical_determinant(w)-rad)==0
    cap = 16*(X+7)
    v = s.Symbol('v')
    for b in K:
        m = adjacent(w,b)*adjacent(g,b).T/cap
        expected = v**2*(v-(X+2)/(X+7))**2 if b<2 else v**3*(v-10/(X+7))
        assert s.factor(m.charpoly(v).as_expr()-expected)==0
    # Fourier transform each high-bit pair in its low bit, then a diagonal order-eight scaling.
    phase = ((-7,-7,-6,-10),(3,3,4,0),(3,3,4,0),(3,3,4,0))
    for key in product(K,repeat=4):
        val = s.factor(sum((-1)**sum((i%2)*j for i,j in zip(key,bits))*
                         w.get(tuple(2*(i//2)+j for i,j in zip(key,bits)),0)
                         for bits in product(range(2),repeat=4)))
        if val:
            exp = sum(phase[h][i] for h,i in enumerate(key))
            assert exp%8==(2 if p in val.free_symbols else 0)


def main():
    a,b,c,d = s.symbols('a b c d')
    # Keep even the zero entries after specialization when differentiating.
    w = universal(a,b,c,d)
    g = alternating_gradient(w)
    cap = sum(w[k]*g[k] for k in w)/4
    assert s.factor(cap-16*(a**2*b**2+a**2*d**2+2*a*b*c*d+b**2*c**2+7*c**2*d**2))==0
    z = s.symbols('z:4')
    for i,vals in enumerate([(1,b,0,1),(a,-a,1,1)]):
        sub = dict(zip((a,b,c,d),vals))
        W,G = [{k:s.factor(f.xreplace(sub)) for k,f in tensor.items()} for tensor in (w,g)]
        if i==0:
            first_curve(b,z,W,G)
        else:
            second_curve(a,z,W,G)
    print('PASS: both kernel pencils, all rank-twelve isomorphisms including null-invariant fibers')


if __name__=='__main__':
    main()
