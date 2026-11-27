"""Solve the ten singly missing supports over characteristic zero, including finite Smith factors."""
from itertools import combinations, product
from functools import reduce
from operator import xor

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from verify_klein_mixed_supports import smith, character
from verify_klein_single_missing_supports import K, SURVIVORS, decode, rank


def blocks(w,p):
    for color in K:
        def entry(i,k):
            key = (i,i^color,k,k^color)
            key = key[-p:]+key[:-p] if p else key
            return w.get(key,0)
        yield color,sp.Matrix([[entry(i,k) for k in K] for i in K])


def equations(keys,variables):
    w = dict(zip(keys,variables))
    result = set()
    for p in (0,1):
        for color,mat in blocks(w,p):
            r = rank(p,color)
            if not r:
                continue
            for i,j in product(combinations(K,r+1),repeat=2):
                value = mat[i,j].det()
                if value==0:
                    continue
                factors = [f for f,k in sp.factor_list(value)[1] if not f.is_Symbol]
                assert len(factors)==1
                result.add(factors[0])
    return sorted(result,key=str)


def normal_forms(mask):
    keys = sorted(decode(mask))
    variables = sp.symbols('a:'+str(len(keys)))
    rows,targets,other = [],[],[]
    for f in equations(keys,variables):
        terms = sp.Poly(f,*variables).terms()
        if len(terms)!=2:
            other.append(f)
            continue
        (a,u),(b,v) = terms
        rows.append([i-j for i,j in zip(a,b)])
        targets.append(-v/u)
    matrix = sp.Matrix(rows)
    d,l,r = smith(matrix)
    rk = matrix.rank()
    assert all(character(row,targets)==1 for row in l.tolist()[rk:])
    gauge = sp.Matrix([[int(key[p]==i) for p,i in product(K,repeat=2)] for key in keys])
    assert matrix*gauge==sp.zeros(matrix.rows,16)
    g = (r.inv()*gauge)[rk:,:]
    _,gl,_ = smith(g)
    gr = g.rank()
    t = sp.symbols('T:'+str(g.rows-gr))
    tail = [character(row,[sp.Integer(1)]*gr+list(t)) for row in gl.inv().tolist()]
    roots = []
    for i,row in enumerate(l.tolist()[:rk]):
        n = int(d[i,i])
        assert n in (1,2,4)
        first = sp.root(character(row,targets),n)
        factors = (1,) if n==1 else (1,-1) if n==2 else (1,-1,sp.I,-sp.I)
        roots.append([first*z for z in factors])
    for choice in product(*roots):
        values = [character(row,list(choice)+tail) for row in r.tolist()]
        sub = dict(zip(variables,values))
        eq = set()
        for f in other:
            val = sp.factor(f.xreplace(sub))
            if val==0:
                continue
            factors = [z for z,k in sp.factor_list(sp.fraction(val)[0])[1] if not z.is_Symbol]
            if not factors:
                eq.add(sp.Integer(1))
            else:
                assert len(factors)==1
                eq.add(factors[0])
        yield dict(zip(keys,values)),t,sorted(eq,key=str)


def chart(n,t):
    v = list(t)
    if n==0:
        v[1]=v[0]*v[3]**2+v[3]
        v[2]=v[0]*v[3]/(v[0]-v[1])
    elif n==2:
        v[2]=1+1/v[1]
    elif n==3:
        v[0]=1/v[1]-1/v[2]
    elif n==4:
        v[0]=v[5]/(v[3]*(v[5]-1))
        v[4]=1-v[0]*v[3]**2/(v[2]*v[5])
        v[1]=1+v[2]*v[4]/v[0]
    elif n==5:
        v[2]=1+v[0]*v[1]
    elif n==6:
        v[0]=v[4]*(v[1]*v[2]-1)/(v[1]-1)
        v[3]=(1-v[4])/(v[2]-v[4])
    elif n==7:
        v[2]=v[3]*v[4]*v[5]-v[3]+1
        v[0]=(v[2]-v[5])/(v[4]-1)
    elif n==8:
        v[4]=v[1]*v[5]+v[2]*v[5]-v[1]*v[2]*v[5]**2
        v[6]=(v[1]*v[3]**2-v[4])/(v[1]*v[3]**2-1)
        v[0]=v[1]*v[3]*(v[2]-v[6])
    elif n==9:
        v[0]=-v[3]**2
        v[1]=1/(v[0]-1)
        v[2]=-2*v[0]**3*v[1]/v[3]**3
    return dict(zip(t,v))


def relations(w,p):
    rel = []
    for color,m in blocks(w,p):
        assert m.rank()==rank(p,color)
        for col in m.columnspace():
            rel.append({(i,i^color):v for i,v in enumerate(col) if v})
    return rel


def degree_four_zero(w):
    cols = [key for key in product(K,repeat=4) if reduce(xor,key)==0]
    rows = []
    for p in range(3):
        for r in relations(w,p):
            for a in product(K,repeat=p):
                for b in product(K,repeat=2-p):
                    if reduce(xor,a+b+next(iter(r))):
                        continue
                    data = {a+k+b:v for k,v in r.items()}
                    rows.append([data.get(key,0) for key in cols])
    return len(cols)-DomainMatrix.from_Matrix(sp.Matrix(rows)).to_field().rank()


def regular_window(t):
    l = 1+t**2
    table = [
        [(-2*t**2/l,-2*t**2,1,1),(t**2,1,-t**2,-1/l),(1,-t**2,-t**2/l,1)],
        [(-2*t**3/l,-2*t,t,1/t),(-1,1,-t**2,t**2/l),(1,1,-t**4/l,-1)],
        [(0,0,0,0),(0,0,-t**2*l,1),(0,0,-t**2,l)],
        [(0,0,0,0),(1,1,0,0),(1,1,0,0)],
    ]
    return {(j^b,j,j^e,j^e^b):sp.sympify(v) for e,group in enumerate(table)
            for b,values in enumerate(group,1) for j,v in enumerate(values) if v}


def pbw(w,u):
    def partition(p):
        return ((0,1),(2,3)) if p%4<2 else ((2,3),(0,1))
    def unit(v):
        num,den = sp.fraction(sp.factor(v))
        return num!=0 and all(f in (u,u**2+1) for z in (num,den) for f,_ in sp.factor_list(z)[1])
    for p in K:
        a,b = partition(p)
        c,d = partition(p+1)
        mixed = []
        counts = [0,0]
        for r in relations(w,p):
            if set(r)<=set(product(a,c)):
                counts[0]+=1
                assert unit(sp.Matrix([[r.get((i,j),0) for j in c] for i in a]).det())
            elif set(r)<=set(product(b,d)):
                counts[1]+=1
                assert unit(sp.Matrix([[r.get((i,j),0) for j in d] for i in b]).det())
            else:
                assert set(r)<=set(product(a,d))|set(product(b,c))
                mixed.append([r.get(key,0) for key in product(b,c)])
        assert counts==[1,1] and unit(sp.Matrix(mixed).det())


def main():
    inv = sp.Symbol('inv')
    for n,mask in enumerate(SURVIVORS):
        live = 0
        for w,t,eq in normal_forms(mask):
            sub = chart(n,t)
            if any(sp.factor(f.xreplace(sub))!=0 for f in eq):
                # All other torsion choices are empty on the coefficient torus.
                assert list(sp.groebner(eq+[inv*sp.prod(t)-1],inv,*t,extension=sp.I))==[1]
                continue
            live += 1
            # On the rational function field in the free coordinates, the only solution is the stated graph.
            dependent = [q for q in t if sub[q]!=q]
            independent = [q for q in t if q not in dependent]
            field = sp.QQ.frac_field(*independent)
            if dependent:
                gb = sp.groebner(eq+[inv*sp.prod(dependent)-1],inv,*dependent,domain=field)
                for q in dependent:
                    assert gb.reduce(sp.cancel(q-sub[q]))[1]==0
            normal = {k:sp.factor(v.xreplace(sub)) for k,v in w.items()}
            value = degree_four_zero(normal)
            assert value==(13 if n==5 else 11 if n==9 else 12)
            if n==9:
                table = regular_window(t[3])
                assert table.keys()==normal.keys() and all(sp.factor(v-normal[k])==0 for k,v in table.items())
                pbw(normal,t[3])
        assert live==1
    print('PASS: ten Laurent rank systems, all Smith root choices, exact quartic character dimensions and the PBW family')


if __name__=='__main__':
    main()
