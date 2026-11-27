"""Exact evaluation-defect ideals at z=-1,0 and the z=0 identification with rectangle 112."""
from itertools import combinations,product
import os
import subprocess
import sympy as s

from verify_exceptional_exchange_sheets import data
from verify_negative_component_quotient import K,a,b,c,d,z,coefficients
from verify_parallel_invertible_rectangles import completed


def ideal(fs):
    return 'ideal('+','.join(str(f).replace('**','^') for f in fs)+')'


def run(z0):
    A = s.Matrix([[a,b],[c,d]])
    D = A.det()
    L = -A.inv().T*s.Matrix([[1,z0],[1,1]])
    N = -A.inv()*s.Matrix([[1,1],[1,z0]])
    y = s.symbols('y:4')
    r = s.Symbol('r')
    args,lines = data()
    assert args==(a,b,c,d,z)
    script = ['ring R=0,(a,b,c,d,r,u),dp;short=0;','ideal I,J;']
    unit = s.Symbol('u')*r*D-1
    for j,(right,low,high) in enumerate(product(range(2),repeat=3)):
        axis = [0]*4
        axis[low],axis[high+2] = 1,r
        (i,k,m,n),(p,q,g,h) = (y,axis) if right else (axis,y)
        pol = [*(D*(s.Matrix([i*p,k*q])-L*s.Matrix([m*g,n*h]))),
               *(D*(s.Matrix([i*q,k*p])-N*s.Matrix([m*h,n*g]))),
               i*g+k*h+m*p+z0*n*q,i*h+k*g+m*q+n*p]
        E = s.Matrix(pol).jacobian(y).applyfunc(s.cancel)
        minors = {s.factor(E.extract(rows,cols).det()) for rows in combinations(range(6),3)
                  for cols in combinations(K,3)}
        sub,t = lines[j]
        expected = [(v-f).subs(z,z0) for v,f in sub.items()]+[(r*r-t).subs(z,z0)]
        if z0==0 and j not in (0,2,4,5):
            expected = [s.S.One]
        script += [f'I=std({ideal([unit]+list(minors))});',f'J=std({ideal([unit]+expected)});',
                   'print(size(reduce(I,J)));print(size(reduce(J,I)));']
    out = subprocess.run([os.environ.get('SINGULAR_BIN','Singular'),'-q'],
                         input='\n'.join(script+['quit;']),text=True,capture_output=True,check=True)
    lines = [l for l in out.stdout.splitlines() if l and not l.startswith('//')]
    assert lines==['0']*16,lines


def rectangular_limit():
    W = coefficients(s.Matrix([[a,b],[c,d]]),s.S.Zero)
    target = completed(112,-1,1,1,a,c,-d,-b)
    sign = ((1,0,0,1),(0,1,0,0),(0,0,1,0),(1,0,0,0))
    for k,f in W.items():
        q = k[0],k[1]^1,k[2],k[3]^1
        assert s.factor(f*(-1)**sum(v[i] for v,i in zip(sign,q))-target.get(q,0))==0


def minus_one():
    W = coefficients(s.Matrix([[a,b],[c,d]]),s.S.NegativeOne)
    for h in (2,3):
        assert s.Matrix(4,4,lambda i,j:W[i,j,i^h,j^h]).rank()==2
    A = s.Matrix([[a,b],[c,d]])
    from verify_klein_elliptic_completion import gram_twist
    m = gram_twist(A,A,s.S.One,s.S.One)[0]
    P = s.Matrix(1,4,[0,0,1,1])
    assert s.Matrix.vstack(*[P*s.Matrix(4,4,m[:,j]) for j in range(16)]).rank()==2


if __name__=='__main__':
    run(-1)
    run(0)
    rectangular_limit()
    minus_one()
    print('PASS: all sixteen localized defect ideals, exact rectangular limit and z=1 obstruction')
