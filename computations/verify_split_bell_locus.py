"""Exact primary-component check for all semisimple invertible alternating tensors."""
from functools import lru_cache
from itertools import permutations, product
from subprocess import run

import sympy as sp

from verify_split_inverse_windows import BITS, bell, bell_coordinates, cut


GROUP = tuple((p,e) for p in permutations(range(4)) for e in product((1,-1),repeat=4)
              if e[0] == 1 and sp.prod(e) == 1)
IDENTITY = (tuple(range(4)), (1,1,1,1))


def multiply(g,h):
    p,e = g
    q,f = h
    signs = tuple(e[i]*f[p[i]] for i in range(4))
    return tuple(q[i] for i in p), tuple(signs[0]*v for v in signs)


@lru_cache(None)
def inverse(g):
    return next(h for h in GROUP if multiply(g,h) == IDENTITY)


def classes():
    remaining = set(GROUP)
    while remaining:
        g = min(remaining)
        orbit = {multiply(h,multiply(g,inverse(h))) for h in GROUP}
        remaining -= orbit
        yield g,len(orbit)


def verify():
    a,b,c,d = x = sp.symbols('a b c d')
    f = -16*cut(bell(x),1).det()
    grad = [f.diff(i) for i in x]
    nonzero = [b,c,d,a-b-c-d,a-b+c+d,a+b-c+d,a+b+c-d]
    quadrics = (a*b+c*d)*(a*c+b*d)*(a*d+b*c)
    lines = [(x[i]-sign*x[j], x[k]-sign*x[l])
             for i,j,k,l in ((0,1,2,3),(0,2,1,3),(0,3,1,2)) for sign in (1,-1)]
    sextic = (d*d+1)**3-4*d**3
    exceptional = [sextic,2*c-d**3-d**2-3*d+1,
                   2*b-d**5+d**4-3*d**3+7*d**2-6*d+2]
    assert sp.Poly(sextic,d).is_irreducible
    def fmt(expr):
        return str(sp.expand(expr.subs(a,1))).replace('**','^')
    def ideal(name,exprs):
        return 'ideal '+name+'='+','.join(map(fmt,exprs))+';'
    header = 'LIB "primdec.lib";ring r=0,(b,c,d),dp;short=0;'
    header += ideal('U',nonzero)+ideal('X',exceptional)+'X=std(X);'
    for n,pair in enumerate(lines):
        header += ideal('L'+str(n),pair)
    header += ideal('surface',[quadrics])
    counts = []
    for index,((p,e),size) in enumerate(classes()):
        equations = [x[i]*e[i]*grad[p[i]]-x[0]*grad[p[0]] for i in (1,2,3)]
        code = header+ideal('I',equations)
        code += ('list parts=minAssGTZ(I);int j,k,ok,good,count,isolated;ideal J;'
                 'for(j=1;j<=size(parts);j++){J=std(parts[j]);ok=1;'
                 'for(k=1;k<=size(U);k++){if(reduce(U[k],J)==0){ok=0;break;}}'
                 'if(ok){count++;good=(size(reduce(surface,J))==0);')
        for n in range(len(lines)):
            code += 'if(size(reduce(L'+str(n)+',J))==0){good=1;}'
        code += ('if(good==0){if(size(reduce(X,J))+size(reduce(J,X))==0)'
                 '{isolated++;}else{print("UNCLASSIFIED");}}}}'
                 'print("COUNT "+string(count)+" "+string(isolated));quit;')
        out = run(['Singular','-q'],input=code,text=True,capture_output=True,check=True).stdout
        assert 'UNCLASSIFIED' not in out and '?' not in out, (index,out)
        line = next(row for row in out.splitlines() if row.startswith('COUNT '))
        counts.append((size,*map(int,line.split()[1:])))
        # Verify dense inclusion in the three stated pieces as well as the upper cover.
        if index == 1:
            assert all(sp.expand(poly.subs({a:1,b:1,c:d})) == 0 for poly in equations)
        if index == 7:
            assert all(sp.expand(poly.subs({a:1,b:-c*d})) == 0 for poly in equations)
    assert counts == [(3,0,0),(1,6,0),(12,0,0),(12,2,0),(32,2,1),
                      (6,1,0),(3,0,0),(3,3,0),(12,0,0),(12,3,0)]
    field_sub = {a:1,b:(d**5-d**4+3*d**3-7*d**2+6*d-2)/2,c:(d**3+d**2+3*d-1)/2}
    def rem(value):
        return sp.rem(sp.expand(value.subs(field_sub)),sextic,d)
    assert rem(quadrics) == -4
    assert rem(sum(i*i for i in x)) != 0
    assert all(rem(x[i]**2-x[j]**2) != 0 for i in range(4) for j in range(i))
    assert not any(all(rem(e[i]*x[i]*x[p[i]]-x[0]*x[p[0]]) == 0 for i in (1,2,3))
                   for p,e in GROUP)
    roots = [sp.sympify(v).subs(field_sub) for v in x]
    reciprocal = bell([sp.invert(v,sextic,d) for v in roots])
    rotated = bell_coordinates({q:reciprocal[q[1:]+q[:1]] for q in BITS})
    assert any(all(sp.rem(roots[p[0]]*rotated[i]-rotated[0]*roots[p[i]]*e[i],sextic,d)==0
                   for i in (1,2,3)) for p,e in GROUP)



if __name__ == '__main__':
    verify()
    print('All 96 comparisons give the three quadrics, six lines and two isolated projective orbits.')
