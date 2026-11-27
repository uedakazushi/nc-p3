"""Repeated surface spectra, monomial isomorphisms, periods and component incidence."""
from collections import defaultdict
from itertools import permutations, product
from subprocess import run

import sympy as sp

from verify_split_spectral_line import F, I, Z, multiply, reconstruct, torus
from verify_split_two_eigenvalue_jordan import bell, character_kernel, decoration, labeled

S, C = sp.symbols('s c', nonzero=True)


def surface(s, c, flip=True, nil=0):
    p = bell((sp.Integer(1), -s*s, s, s))
    p[0,1,1,0] = nil
    maps = multiply(torus(c, 1/c), multiply(F if flip else (I,)*4, (Z,I,Z,I)))
    return labeled(reconstruct(p,maps))


def split_coordinates(w):
    def f(q):
        return w[tuple(map(int,q))]
    return (sp.factor(f('0213')*f('1302')/(f('0312')*f('1203'))),
            sp.factor(-f('2031')*f('0132')*f('1230')/(f('0231')*f('2130')*f('1032'))))


def exponents(value):
    value = sp.factor(value)
    num, den = sp.fraction(value)
    constant = sp.Poly(num,S,C).LC()/sp.Poly(den,S,C).LC()
    powers = sp.factor(value/constant).as_powers_dict()
    factors = (S, S-1, S+1, S**2+1, C)
    exps = [powers.get(f,0) for f in factors]
    two = next(j for j in range(-20,21) if abs(constant)==sp.Rational(2)**j)
    assert sp.cancel(value-constant*sp.prod(f**p for f,p in zip(factors,exps))) == 0
    return exps+[two,int(bool(constant<0))]


def diagonal_equivalence(ker, ratios):
    residual = ker.T*sp.Matrix([exponents(r) for r in ratios])
    assert residual[:,:-1].is_zero_matrix and all(v%2==0 for v in residual[:,-1])


def isomorphisms(w):
    words = tuple(sorted(w))
    ker = character_kernel(words)
    f = (S**2-1)/(2*S)
    assert all(sp.cancel(x-y)==0 for x,y in zip(split_coordinates(w),(f**2,f*C)))
    from sympy.matrices.normalforms import hermite_normal_form
    lattice = ker.T*sp.Matrix([exponents(w[q]) for q in words])[:,:5]
    assert hermite_normal_form(lattice.T).T.tolist()==[[-2,2,2,0,0],[-2,0,0,2,0],[-1,1,1,0,1]]
    for shift in (0,1):
        pairs = set()
        for p in permutations(range(4)):
            def target(q):
                r = tuple(p[j] for j in q)
                return r[shift:]+r[:shift]
            if any(target(q) not in w for q in words):
                continue
            new = {q:w[target(q)] for q in words}
            u,v = split_coordinates(new)
            t = next(x for x in (S,(S+1)/(S-1)) if sp.cancel(((x*x-1)/(2*x))**2-u)==0)
            d = sp.factor(2*t*v/(t*t-1))
            pairs.add((t==S, d))
            ref = {q:sp.factor(a.subs({S:t,C:d},simultaneous=True)) for q,a in w.items()}
            diagonal_equivalence(ker,[sp.factor(new[q]/ref[q]) for q in words])
        sign = (-1)**shift
        assert pairs==set(product((True,False),(sign*C,sign/C)))
    for sign in (sp.Integer(1),sp.Integer(-1)):
        decorated = surface(S,sign,nil=1)
        assert decoration(decorated,(0,3)) == -1
        assert decoration(decorated,(3,0)) == -1
        rotated = {q[2:]+q[:2]:v for q,v in decorated.items()}
        words = tuple(sorted(decorated))
        diagonal_equivalence(character_kernel(words),[sp.factor(sp.sympify(rotated[q])/decorated[q]) for q in words])


def tangent_blocks(w):
    pairs = sorted(product(range(4),repeat=2),key=lambda q:q[0]<=q[1])
    kernels = []
    for i in range(4):
        def val(p):
            return w.get(p[-i:]+p[:-i] if i else p,0)
        matrix, pivots = sp.Matrix([[val(p+b) for p in pairs] for b in pairs]).rref()
        assert pivots == tuple(range(6))
        kernels.append([{q:(-matrix[j,k] if j<6 else 1) for j,q in enumerate(pairs)
                         if (matrix[j,k] if j<6 else j==k)} for k in range(6,16)])
    def character(q):
        return q.count(0)-q.count(2),q.count(1)-q.count(3)
    words, blocks = defaultdict(list),defaultdict(list)
    for q in product(range(4),repeat=4):
        words[character(q)].append(q)
    for i in (0,1):
        for a,b in product(kernels[i],kernels[i+2]):
            row = {}
            for x,u in a.items():
                for y,v in b.items():
                    q = x+y
                    row[q[-i:]+q[:-i] if i else q] = sp.factor(u*v)
            h = character(next(iter(row)))
            assert all(character(q)==h for q in row)
            den = sp.lcm([sp.denom(v) for v in row.values()])
            assert all(sp.rem((S*C*(S**4-1))**10, f, S, C)==0
                       for f,_ in sp.factor_list(den)[1])
            blocks[h].append([str(sp.factor(row.get(q,0)*den)).replace('**','^') for q in words[h]])
    return [(len(words[h]),rows) for h,rows in blocks.items()]


def check_tangents(w):
    blocks = tangent_blocks(w)
    script = ['LIB "jacobson.lib";ring r=(0,c),s,dp;minpoly=MODULUS;poly support=s*(s^4-1);',
              'proc clean(poly f){poly d=gcd(f,support);while(deg(d)>0){f=f/d;d=gcd(f,support);}return(f/leadcoef(f));}']
    for size, rows in blocks:
        values = ','.join(v for row in rows for v in row)
        script.append(f'matrix m[{len(rows)}][{size}]={values};')
        if size==1:
            script.append('poly d=m[1,1];int j;for(j=2;j<=nrows(m);j++){d=gcd(d,m[j,1]);};'
                          'if(d==0 || deg(clean(d))>0){print("JUMP");};kill d,j,m;')
        else:
            script.append(f'matrix h=smith(transpose(m));int j;int z=0;for(j=1;j<={size};j++)'
                          '{if(j>ncols(h)){z++;}else{if(h[j,j]==0){z++;}else{if(deg(clean(h[j,j]))>0)'
                          '{print("JUMP");};};};};print("ZERO "+string(z));kill j,z,h,m;')
    for modulus in ('c^2+1','c^2+c+1','c^2-c+1','c^4+1'):
        out = run(['Singular','-q'],input='\n'.join(script+['quit;']).replace('MODULUS',modulus),
                  text=True,capture_output=True,check=True).stdout
        assert '?' not in out and 'JUMP' not in out, out
        assert sum(int(z.split()[1]) for z in out.splitlines() if z.startswith('ZERO')) == 65


if __name__=='__main__':
    tensor = surface(S,C)
    isomorphisms(tensor)
    check_tangents(tensor)
    print('Repeated surface classes, all shifts and all required tangent ranks verified.')
