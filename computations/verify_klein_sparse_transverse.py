"""Solve all sparse transverse doublet supports on their coefficient tori."""
from itertools import combinations, permutations, product
from pathlib import Path
from subprocess import check_output
from tempfile import TemporaryDirectory

import sympy as s

from verify_klein_three_lines import K, WORDS, entry, smith, twist
from verify_klein_mixed_positive import bit, fourth_character
from verify_klein_two_doublet_supports import candidates

POINT = 0x0ff066663c3caa55
CURVE = 0x0ff066667e7eaff5


def torus(mask):
    words = [w for w in WORDS if mask&bit(w)]
    index = {w:i for i,w in enumerate(words)}
    equations,others = {},{}
    for p,b in product(range(2),K):
        size = 3 if b in (0,1 if p==0 else 2) else 2
        for I,J in product(combinations(K,size),repeat=2):
            terms = []
            for pi in permutations(range(size)):
                keys = [entry(p,b,i,J[j]) for i,j in zip(I,pi)]
                if all(k in index for k in keys):
                    sign = (-1)**sum(i>j for i,j in combinations(pi,2))
                    terms.append((tuple(sorted(index[k] for k in keys)),sign))
            if not terms:
                continue
            assert len(terms)>1
            common = set.intersection(*(set(key) for key,_ in terms))
            terms = sorted((tuple(i for i in key if i not in common),v) for key,v in terms)
            if len(terms)==2:
                (k,a),(l,beta) = terms
                row = tuple(int(i in k)-int(i in l) for i in range(len(words)))
                assert row not in equations or equations[row]==-a*beta
                equations[row]=-a*beta
            else:
                others[tuple((key,v*terms[0][1]) for key,v in terms)]=True
    m = s.Matrix(list(equations))
    d,left,right = smith(m)
    ds = [n for n in d.diagonal() if n]
    assert all(abs(n)==1 for n in ds)
    rk = len(ds)
    signs = left*s.Matrix([int(v==-1) for v in equations.values()])
    assert all(e%2==0 for e in signs[rk:])
    gauge = s.Matrix([[int(w[p]==i) for p,i in product(K,repeat=2)] for w in words])
    assert m*gauge==s.zeros(m.rows,16)
    e,rows,_ = smith((right.inv()*gauge)[rk:,:])
    gr = sum(bool(a) for a in e.diagonal())
    exponents = right[:,rk:]*rows.inv()[:,gr:]
    y = s.symbols('y:'+str(exponents.cols))
    sign_vec = right[:,:rk]*signs[:rk,:]
    w = {key:s.Integer(-1)**int(sign_vec[j]%2)*s.prod(v**int(e) for v,e in zip(y,row))
         for j,(key,row) in enumerate(zip(words,exponents.tolist()))}
    minors = set()
    for f in others:
        poly = s.factor(sum(v*s.prod(w[words[j]] for j in key) for key,v in f))
        if poly==0:
            continue
        if not y:
            return w,y,[s.Integer(1)]
        unit_free = s.Poly(s.fraction(poly)[0],*y).terms_gcd()[1].monic()
        if len(unit_free.terms())==1:
            return w,y,[s.Integer(1)]
        minors.add(unit_free.as_expr())
    return w,y,sorted(minors,key=str)


def primes(eq,y):
    if eq==[1]:
        return []
    if not eq:
        return [[]]
    script = ('LIB "primdec.lib"; LIB "elim.lib"; ring r=0,('+','.join(map(str,y))+'),dp; ideal J='+
              ','.join(str(f).replace('**','^') for f in eq)+'; list S=sat(J,ideal('+
              '*'.join(map(str,y))+')); ideal B=S[1]; if(reduce(1,std(B))!=0){list L=minAssGTZ(B);'+
              'int j; for(j=1;j<=size(L);j++){print("PRIME");print(L[j]);};}; print("END"); quit;')
    with TemporaryDirectory() as temp:
        path = Path(temp)/'torus.sing'
        path.write_text(script)
        output = check_output(['Singular','-q',str(path)],text=True)
    assert 'END' in output and '?' not in output
    return [[s.sympify(line.rstrip(',').replace('^','**')) for line in block.splitlines() if line]
            for block in output.split('END')[0].split('PRIME\n')[1:]]


def rational_chart(prime):
    pool,sub = list(prime),{}
    while pool:
        choices = []
        for f in pool:
            for x in f.free_symbols:
                poly = s.Poly(f,x)
                if poly.degree()==1:
                    val = s.factor(-poly.nth(0)/poly.nth(1))
                    choices.append((s.count_ops(val),str(x),str(val),x,val))
        _,_,_,x,val = min(choices,key=lambda item:item[:3])
        assert val!=0
        sub = {j:s.factor(a.xreplace({x:val})) for j,a in sub.items()}
        sub[x] = val
        pool = list({s.fraction(s.factor(f.xreplace({x:val})))[0] for f in pool}-{0})
    assert all(s.factor(f.xreplace(sub))==0 for f in prime)
    return sub


def boundary_model(mask,w):
    t = s.Symbol('t',nonzero=True)
    y = s.Symbol('y0')
    gauge = ((1,-1,0,0),(-1,1,0,0),(1,1,0,0),(1,1,0,0)) if mask==CURVE else ((0,)*4,)*4
    w = {key:s.factor(val.subs(y,t**2)*t**sum(gauge[p][key[p]] for p in K)) for key,val in w.items()}
    assert all(val==w[tuple(i^3 for i in key)] for key,val in w.items())
    z = t**-2 if mask==CURVE else s.Integer(1)
    diag = ((1,z,1,-z),(1,z,1,z),(1,1,-1,1))
    S = {}
    for key in WORDS:
        f = sum((-1)**sum((i//2)*b for i,b in zip(key,bits))*
                w.get(tuple((i%2)^(3*b) for i,b in zip(key,bits)),0) for bits in product(range(2),repeat=4))
        a,b,c,d = key
        S[a,b,c,d^1] = s.factor(diag[0][b]*diag[1][c]*diag[2][d^1]*twist(key)*f/2)
    assert all(val==S[perm] for key,val in S.items() for perm in permutations(key))
    P,Q = ((3*t**2+1)/t**4,(1-t**2)/t**4) if mask==CURVE else (3,-1)
    m = [s.Matrix([[S.get((i,i^b,k,k^b^1),0) for k in K] for i in K]) for b in K]
    a,c = s.Matrix([1,P,1,Q]),s.Matrix([1,Q,1,P])
    expected = [s.Matrix.hstack(a,a,c,c),s.Matrix.vstack(a.T,a.T,c.T,c.T),
                s.Matrix([1,Q,1,Q])*s.ones(1,4),s.ones(4,1)*s.Matrix([[1,Q,1,Q]])]
    assert all((left-right).applyfunc(s.factor)==s.zeros(4) for left,right in zip(m,expected))
    # Each opposite block is full except for y = 1,-1 on the one-dimensional torus.
    for b in K:
        A = s.Matrix([[w.get((i,j,i^b,j^b),0) for j in K] for i in K])
        if mask==POINT:
            assert A.rank()==4
        else:
            det = s.factor(A.det())
            assert det!=0 and set(f for f,_ in s.factor_list(s.fraction(det)[0])[1] if f!=t)=={t-1,t+1,t**2+1}
            assert all(A.subs(t,j).rank()==3 for j in (1,-1,s.I,-s.I))


def main():
    number,good = 0,[]
    for mask in candidates(2):
        w,y,eq = torus(mask)
        for prime in primes(eq,y):
            number += 1
            sub = rational_chart(prime)
            candidate = {key:s.factor(val.xreplace(sub)) for key,val in w.items()}
            if any(fourth_character(candidate,p,0)>9 for p in K):
                continue
            assert not prime and mask in (POINT,CURVE)
            boundary_model(mask,candidate)
            good.append(mask)
    assert number==63 and good==[POINT,CURVE]
    print('PASS: 51 supports, 63 coefficient components, two regular boundary forms, every opposite rank')


if __name__=='__main__':
    main()
