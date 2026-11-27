"""Unrestricted smoothness and Klein-equivariant smoothings of the three rectangular families."""
from functools import reduce
from itertools import combinations, product
from operator import xor
import os
import random
import subprocess

from flint import nmod_mat
from verify_parallel_invertible_rectangles import K, completed, matrix, s
from verify_parallel_horizontal_moduli import window
from verify_parallel_third_rectangle import third_rectangle

PRIME = 1009


def residue(x):
    n,d = s.fraction(x)
    return int(n)*pow(int(d),-1,PRIME)%PRIME


def examples():
    u,v,z,w = map(s.Integer,(2,3,5,4))
    t = (u+v)/(1+u*v)
    a,y = u*(z+v)/(1+u*v),u*(w+v)/(1+u*v)
    return [completed(76,t,u*t,v*t,a,1,1,v*(1-u*z)/(u+v)),
            window(t,u*t,v*t,a,y),third_rectangle(s.Integer(2),s.Integer(3),s.Integer(5))]


def word(p,i,j,h,g):
    key = (i,i^h,j,j^g)
    return key[-p:]+key[:-p] if p else key


def schur_differentials(tensor,dimensions,check_smoothing=True):
    W = {k:residue(v) for k,v in tensor.items()}
    words = [[k for k in product(K,repeat=4) if reduce(xor,k)==h] for h in K]
    eqs = [[] for _ in K]
    for p in (0,1):
        kernels = {}
        for h in K:
            M = nmod_mat(4,4,[W.get(word(p,i,j,h,h),0) for i,j in product(K,repeat=2)],PRIME)
            L,ld = M.transpose().nullspace()
            R,rd = M.nullspace()
            kernels[h] = (L,ld,R,rd)
        for h,g in product(K,repeat=2):
            L,ld,_,_ = kernels[h]
            _,_,R,rd = kernels[g]
            for i,j in product(range(ld),range(rd)):
                row = {word(p,k,l,h,g):int(L[k,i])*int(R[l,j])%PRIME for k,l in product(K,repeat=2)}
                eqs[h^g].append([row.get(k,0) for k in words[h^g]])
    equations = [nmod_mat(rows,PRIME) for rows in eqs]
    assert [64-m.rank() for m in equations]==list(dimensions)
    for h,m in enumerate(equations):
        orbit = nmod_mat([[W.get(k[:p]+(k[p]^h,)+k[p+1:],0) if k[p]==i else 0 for k in words[h]]
                          for p,i in product(K,repeat=2)],PRIME)
        assert orbit.rank()==(13,16,16,16)[h] and (m*orbit.transpose()).rank()==0
    if not check_smoothing:
        return
    basis,rank = equations[0].nullspace()
    rng = random.Random(456)
    coeff = [rng.randrange(PRIME) for _ in range(rank)]
    tangent = {k:sum(c*int(a) for c,a in zip(coeff,row))%PRIME for k,row in zip(words[0],basis.tolist())}
    for p in (0,1):
        for h in K:
            def entry(i,j):
                k = word(p,i,j,h,h)
                return W.get(k,0),tangent[k]
            if h>=2:
                assert all(any(entry(i,j)!=(0,0) for j in K) for i in K)
                assert all(any(entry(i,j)!=(0,0) for i in K) for j in K)
            else:
                for a,b in ((0,2),(2,0)):
                    (v,dv),(w,dw),(x,dx),(y,dy) = [entry(a+i,b+j) for i,j in product(range(2),repeat=2)]
                    assert (v*y-w*x)%PRIME or (dv*y+v*dy-dw*x-w*dx)%PRIME
                    if (p,h,a,b)==(1,0,0,2):
                        assert (v*y+w*x)%PRIME


def point_algebra(rel,hyperplane):
    coords = s.symbols('a b c d')
    M = s.Matrix([[col[k^h]*coords[k^h] for k in K] for h,col in rel])
    def fmt(x):
        return str(x).replace('**','^')
    polys = []
    for I in combinations(range(6),4):
        denom,f = s.Poly(M.extract(I,K).det(method='domain-ge'),*coords).clear_denoms(convert=True)
        assert denom%PRIME
        polys.append(fmt(f.as_expr()))
    form = fmt(sum(i*j for i,j in zip(coords,hyperplane)))
    code = [f'ring r={PRIME},(a,b,c,d),dp;short=0;ideal F='+','.join(polys)+';poly l='+form+';ideal H;']
    for i in K:
        chart = ','.join([str(x) for x in coords[:i]]+[str(coords[i])+'-1'])
        code.append('H=std(F+ideal(l,'+chart+'));print(H[1]);')
    code.extend(['H=std(F+ideal(l-1));print(vdim(H));ideal B=kbase(H);int j;',
                 'for(j=1;j<=size(B);j++){print(string(B[j]));}'])
    for x in coords:
        code.append('for(j=1;j<=size(B);j++){print(string(reduce('+str(x)+'*B[j],H)));}')
    res = subprocess.run([os.environ.get('SINGULAR_BIN','Singular'),'-q'],input='\n'.join(code+['quit;']),
                         text=True,capture_output=True,check=True)
    data = [r for r in res.stdout.splitlines() if r and not r.startswith('//')]
    assert data[:5]==['1']*4+['20'] and len(data)==105
    def poly(line):
        return s.Poly(s.sympify(line.replace('^','**')),*coords)
    ids = {poly(f).monoms()[0]:j for j,f in enumerate(data[5:25])}
    matrices = []
    for i in K:
        rows = []
        for line in data[25+20*i:45+20*i]:
            row = [0]*20
            for m,c in poly(line).terms():
                if c:
                    row[ids[m]]=int(c)%PRIME
            rows.append(row)
        matrices.append(nmod_mat(rows,PRIME).transpose())
    return matrices,ids[(0,0,0,0)]


def conormals(tensor):
    rel = [[(h,col) for h in K for col in matrix(tensor,p,h).columnspace()] for p in K]
    rng = random.Random(772)
    a,unit = point_algebra(rel[0],[rng.randrange(1,30) for _ in K])
    ident = nmod_mat([[int(i==j) for j in range(20)] for i in range(20)],PRIME)
    zero = ident*0
    def partner(p,point,reverse=False):
        eq = [[point[k^h]*residue(col[k if reverse else k^h]) for k in K] for h,col in rel[p]]
        for _ in range(12):
            rows = []
            for _ in range(3):
                weights = [rng.randrange(PRIME) for _ in eq]
                rows.append([sum((m*c for m,c in zip(column,weights)),zero) for column in zip(*eq)])
            rows.append([ident*rng.randrange(1,PRIME) for _ in K])
            M = nmod_mat([[int(rows[i][j][r,t]) for j in K for t in range(20)]
                          for i in K for r in range(20)],PRIME)
            if M.rank()==80:
                break
        else:
            raise AssertionError('missing invertible evaluation chart')
        inv = M.inv()
        other = [nmod_mat([[int(inv[20*j+r,60+t]) for t in range(20)] for r in range(20)],PRIME) for j in K]
        assert all(sum((x*y for x,y in zip(row,other)),zero).rank()==0 for row in eq)
        return other
    def equal(x,y):
        return all((a*b-c*d).rank()==0 for a,c in zip(x,y) for b,d in zip(y,x))
    b = partner(0,a)
    c = partner(1,b)
    d = partner(2,c)
    assert equal(a,partner(3,d)) and equal(d,partner(3,a,True))
    rows = [[int(M[i,unit]) for i in range(20)] for x,y,z,w in product(a,b,c,d) for M in (x*y*z*w,)]
    assert nmod_mat(rows,PRIME).rank()==20


if __name__=='__main__':
    for name,tensor in zip(('76','112','third'),examples()):
        schur_differentials(tensor,(18,18,20,20))
        conormals(tensor)
        print('PASS: unrestricted Schur identities and smoothing of all thin blocks for',name,flush=True)

    for shape in (76,112):
        for sign in (1,-1):
            y = 1 if shape==76 else 4
            vertical = completed(shape,sign,1,1,3*sign,1,y,7*sign)
            dims = (17,18,20,20) if sign==1 else (17,16,16,16)
            schur_differentials(vertical,dims)
            u = v = s.Integer(sign)
            z,w = s.Integer(3),s.Integer(4)
            t = (u+v)/(1+u*v)
            a,y = u*(z+v)/(1+u*v),u*(w+v)/(1+u*v)
            cross = (completed(76,t,1,1,a,1,1,(1-a)/t**2) if shape==76
                     else window(t,1,1,a,y))
            schur_differentials(cross,(19,18,20,20),False)
    print('PASS: both vertical components and both signs of the component intersections')
