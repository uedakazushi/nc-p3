"""Stabilizers, normalized orbits, and rotations of the 76-type rectangle."""
from itertools import product

from verify_parallel_invertible_rectangles import K,T,X,bit,candidates,completed,s,smith,template,top
from verify_third_rectangle_orbits import support_maps


def check():
    t,b,x,d = s.symbols('t b x d',nonzero=True)
    m = t*t
    c = (m-b)/(1-b)
    D = c*(1-x)/m
    w = completed(76,t,b,c,x,1,1,D)
    holes = [k for k,z in template(76).items() if z in X]
    thin = [k for k in w if k not in holes]
    def weight(k):
        return [int(k[p]==j) for p,j in product(K,repeat=2)]
    def exponent(f):
        q = s.factor(t*s.diff(f,t)/f)
        assert q.is_Integer and f==t**q
        return q
    G = s.Matrix([weight(k)+[exponent(w[k])] for k in thin])
    step = [(t,b,c,x),(1/t,c/m,1/b,c*(1-x)/m),(t,m/b,m/c,1-c+c*x/b),(1/t,1/c,b/m,1-x/b)]
    for idx in (85,82):
        words = [k for k in w if candidates(1)[idx]&bit(k)]
        diag = smith(s.Matrix(list(map(weight,words))))[0]
        assert [n for n in diag.diagonal() if n]==[1]*11+[2,2]
        for r in K:
            rotated = [k[r:]+k[:r] for k in words]
            perms = list(support_maps(rotated,words))
            assert len(perms)==(2 if idx==85 else int(r%2==0))
            if idx==85 and r==0:
                assert perms[1]==tuple(tuple(i^h for i in K) for h in (2,3,2,3))
            vertical_results = []
            for pi in perms:
                def move(raw):
                    return {tuple(p[i] for p,i in zip(pi,k[r:]+k[:r])):v for k,v in raw.items()}
                z = move(w)
                sol = next(iter(s.linsolve((G,s.Matrix([exponent(z[k]) for k in thin])))))
                subs = dict.fromkeys(set().union(*(v.free_symbols for v in sol)),0)
                es = [v.xreplace(subs) for v in sol]
                def normalize(raw,tt):
                    source = move(raw)
                    def val(k):
                        return source[k]*tt**-sum(es[4*p+i] for p,i in enumerate(k))
                    assert all(s.factor(val(k)-template(76)[k].subs(T,tt**es[-1]))==0 for k in thin)
                    u,v = val(holes[4]),val(holes[5])
                    return tuple(s.factor(val(k)/(u if j in (0,1,2,4) else v)) for j,k in enumerate(holes[:8]))
                def compare(left,right):
                    assert all(s.factor(f-g)==0 for f,g in zip(left,right))
                if idx==85:
                    tt,bb,cc,xx = step[r]
                    compare((t**es[-1],),(tt,))
                    compare(normalize(w,t),top(76,tt,bb,cc,xx,1,1,cc*(1-xx)/tt**2))
                for t0 in (s.S.One,-s.S.One):
                    xx = x if idx==85 else s.S.One
                    raw = completed(76,t0,1,1,xx,1,1,d)
                    vec = normalize(raw,t0)
                    compare((t0**es[-1],),(t0,))
                    A,DD = vec[2],vec[6]
                    compare(vec,top(76,t0,1,1,A,1,1,DD))
                    if t0==1:
                        vertical_results.append((A,DD))
                    else:
                        compare((A,DD),vertical_results[-1])
            expected = ([(x,d),(1-d,1-x)] if r%2==0 else [(d,x),(1-x,1-d)]) if idx==85 else (
                        [(s.S.One,d)] if r%2==0 else [])
            assert set(vertical_results)==set(expected)


if __name__=='__main__':
    check()
    print('PASS: all support permutations, normalized coefficients, shifts and diagonal stabilizers on 82 and 85')
