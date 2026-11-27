"""Weyl action on all simple surface sheets and tangent ranks on the extra curve."""
from collections import defaultdict
from functools import reduce
from itertools import product
from subprocess import run
import re

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from verify_split_inverse_windows import BITS, bell, transport
from verify_split_spectral_line import F, I, X, Z, compare, gauge, multiply, reconstruct

U,V = sp.symbols('u v',nonzero=True)
PAULI = (I,X,Z,X*Z)
D0 = (Z,I,Z,I)
EXTRA = (I,X,Z,X*Z)


def surface(u,v):
    return bell((sp.Integer(1),-u*v,u,v))


def equal_lines(a,b):
    p = BITS[0]
    return all(sp.cancel(a[q]*b[p]-b[q]*a[p])==0 for q in BITS)


def sheets():
    p = surface(U,V)
    phase = sp.diag(1,sp.I)
    choices = [((-U,V),(I,I,X*Z,X*Z),(3,2,1,0)),
               ((1/U,V),(I,X,X,I),(1,0,3,2)),
               ((V,U),(phase,phase.inv(),phase.inv(),phase),(2,1,0,3))]
    maps = [multiply((h,)*4,D0) for h in PAULI]
    for (x,y),basis,permutation in choices:
        assert equal_lines(transport(p,basis,range(4)),surface(x,y))
        for j,idx in enumerate(permutation):
            compare(gauge(maps[j],basis),maps[idx])
    phi = lambda x:(x-1)/(x+1)
    assert equal_lines({q[1:]+q[:1]:v for q,v in p.items()},surface(-1/phi(V),phi(U)))
    for h,m in enumerate(maps):
        compare(tuple(m[(i-1)%4].inv() for i in range(4)),maps[h^2])
    p = surface(sp.I,V)
    assert equal_lines(transport(p,EXTRA,range(4)),p)
    extra = multiply(EXTRA,D0)
    compare(gauge(extra,EXTRA),multiply((Z,)*4,extra))
    swap = tuple(extra[i^1].T for i in range(4))
    compare(gauge(swap,(I,Z,Z,I)),multiply(F,extra))
    return reconstruct(p,extra)


def tangent_blocks(w):
    field = sp.QQ_I.frac_field(V)
    tensor = {q:field.from_sympy(sp.sympify(v)) for q,v in w.items()}
    labels = [[(j^(i%2)) if j<2 else (j^int(i>=2)) for j in range(4)] for i in range(4)]
    kernels = []
    for i in range(4):
        pairs = sorted(product(range(4),repeat=2),key=lambda q:labels[i][q[0]]<=labels[(i+1)%4][q[1]])
        def val(q):
            return tensor.get(q[-i:]+q[:-i] if i else q,field.zero)
        mat,pivots = DomainMatrix.from_list([[val(p+b) for p in pairs] for b in pairs],field).rref()
        assert pivots == tuple(range(6))
        r = mat.to_list()
        kernels.append([{q:(-r[j][k] if j<6 else field.one) for j,q in enumerate(pairs)
                         if (r[j][k] if j<6 else j==k)} for k in range(6,16)])
    parent = {}
    def root(q):
        while parent.get(q,q)!=q:
            q = parent[q]
        return q
    rows = []
    var = field.gens[0]
    allowed = (var*(var**4-1)).numer
    for i in (0,1):
        for a,b in product(kernels[i],kernels[i+2]):
            row = {}
            for x,u in a.items():
                for y,v in b.items():
                    q = x+y
                    row[q[-i:]+q[:-i] if i else q] = u*v
            den = reduce(lambda a,b:a.lcm(b),(v.denom for v in row.values()))
            assert not allowed**20 % den
            first = root(next(iter(row)))
            for q in row:
                parent[root(q)] = first
            rows.append({q:(v.numer*(den//v.denom)).as_expr() for q,v in row.items()})
    cols, blocks = defaultdict(list),defaultdict(list)
    for q in product(range(4),repeat=4):
        cols[root(q)].append(q)
    for row in rows:
        blocks[root(next(iter(row)))].append(row)
    return [(c,blocks[h]) for h,c in cols.items()]


def verify_tangent(w):
    script = ['LIB "jacobson.lib";ring r=(0,c),v,dp;minpoly=c^2+1;poly support=v*(v^4-1);',
              'proc clean(poly f){poly d=gcd(f,support);while(deg(d)>0){f=f/d;d=gcd(f,support);}return(f/leadcoef(f));}']
    for cols,rows in tangent_blocks(w):
        size = len(cols)
        coeffs = ','.join(str(r.get(q,0)).replace('**','^').replace('I','c') for r in rows for q in cols)
        coeffs = re.sub(r'v\^(\d+)/',r'(v^\1)/',coeffs)
        script.append(f'matrix m[{len(rows)}][{size}]={coeffs};')
        if size==1:
            script.append('poly d=m[1,1];int j;for(j=2;j<=nrows(m);j++){d=gcd(d,m[j,1]);};'
                          'if(d==0 || deg(clean(d))>0){print("JUMP");};kill d,j,m;')
        else:
            script.append(f'matrix h=smith(transpose(m));int j;int z=0;for(j=1;j<={size};j++)'
                          '{if(j>ncols(h)){z++;}else{if(h[j,j]==0){z++;}else{if(deg(clean(h[j,j]))>0)'
                          '{print("JUMP");};};};};print("ZERO "+string(z));kill j,z,h,m;')
    out = run(['Singular','-q'],input='\n'.join(script+['quit;']),text=True,capture_output=True,check=True).stdout
    assert '?' not in out and 'JUMP' not in out, out
    assert sum(int(z.split()[1]) for z in out.splitlines() if z.startswith('ZERO'))==61


if __name__=='__main__':
    verify_tangent(sheets())
    print('All surface sheet generators, exchange, one-shift and extra smooth component verified.')
