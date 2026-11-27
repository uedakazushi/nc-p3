"""All exceptional low-rank evaluations on the four vertical components."""
from itertools import combinations
from pathlib import Path
from subprocess import check_output
from tempfile import TemporaryDirectory

from verify_parallel_invertible_rectangles import A,B,C,D,K,completed,matrix,product,s
from verify_klein_three_lines import twist

SUPPORTS = ((0,2),(0,3),(1,2),(1,3),(0,1,2,3))


def vertical(shape,t):
    w = completed(shape,t,1,1,A,B,C,D)
    untwisted = {k:val*(twist(k) if t==-1 else 1) for k,val in w.items()}
    even = (1,1,-1,1) if t==-1 else (1,1,1,1)
    odd = ((1,-1,1,1) if shape==76 else (-1,1,-1,1)) if t==-1 else (1,1,1,1)
    gram = {(i,j,k,l^1):s.factor(val*even[k]*odd[l]) for (i,j,k,l),val in untwisted.items()}
    assert all(v==gram.get((k,j,i,l),0)==gram.get((i,l,k,j),0) for (i,j,k,l),v in gram.items())
    for h in (2,3):
        m = s.Matrix([[gram.get((i,j,i^h,j^h^1),0) for j in K] for i in K])
        assert not set().union(*(f.free_symbols for f in m)) and m.rank()==2
    windows = []
    for p in (0,1):
        rows = []
        for h in K:
            m = matrix(w,p,h)
            columns = [m[:,0],m[:,1]] if h<2 else m.columnspace()[:1]
            for v in columns:
                rows.append((h,s.Matrix([v[j]*((-1)**((j>>1)*((j^h)&1)) if t==-1 else 1) for j in K])))
        windows.append(rows)
    return windows


def pure_check(w,p):
    thin = []
    for h in (2,3):
        m = matrix(w,p,h)
        thin.append(m.columnspace()[0])
    u,v = thin
    ranks = s.Matrix([[u[0]*u[1],-v[0]*v[1]],[u[2]*u[3],-v[2]*v[3]]])
    assert ranks.det()!=0


def expected(shape,t,p,nz):
    if shape==76:
        if (p,nz)==(0,(0,2)):
            return [D,A-B]
        if (p,nz)==(1,(1,2)):
            return [C-D,A]
    else:
        return {(0,(0,2)):[B-t*D,A-D],(0,(1,2)):[C-t*D,B-t*D],
                (1,(1,2)):[C-t*D,A],(1,(1,3)):[C,A-D]}.get((p,nz),[1])
    return [1]


def projection(window,shape,t,p,nz):
    xs = s.symbols('u:4')
    u = [xs[j] for j in nz[1:]]
    coord = {k:(1 if k==nz[0] else xs[k] if k in nz else 0) for k in K}
    m = s.Matrix([[v[j^h]*coord[j^h] for j in K] for h,v in window])
    eq = {s.factor(m.extract(I,J).det()) for I,J in product(combinations(range(6),3),combinations(K,3))}
    aux = s.Symbol('v')
    polynomials = list(eq-{0})+[aux*s.prod(u)*(A*C-B*D)-1]
    def ideal(polys):
        return ','.join(str(f).replace('**','^') for f in polys)
    ring = ','.join(map(str,[aux,*u,A,B,C,D]))
    eliminate = '*'.join(map(str,[aux,*u]))
    code = (f'LIB "primdec.lib";ring r=0,({ring}),dp;ideal J={ideal(polynomials)};'
            f'ideal L=std(radical(eliminate(J,{eliminate})));ideal H=std(ideal({ideal(expected(shape,t,p,nz))}));'
            'if(size(reduce(L,H))==0 && size(reduce(H,L))==0){print("PASS");};quit;')
    with TemporaryDirectory() as tmp:
        path = Path(tmp)/'projection.sing'
        path.write_text(code)
        out = check_output(['Singular','-q',str(path)],text=True)
    assert 'PASS' in out.splitlines() and '?' not in out, (shape,t,p,nz,out)


if __name__=='__main__':
    for shape,t in product((76,112),(1,-1)):
        windows = vertical(shape,t)
        W = completed(shape,t,1,1,A,B,C,D)
        for p in K:
            pure_check(W,p)
        for p,window in enumerate(windows):
            for nz in SUPPORTS:
                projection(window,shape,t,p,nz)
    print('PASS: both Gram symmetries, all thin determinants and all forty exceptional evaluation projections')
