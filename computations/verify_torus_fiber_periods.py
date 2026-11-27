"""Actions of position shifts on positive-weight slices, including the finite sign stabilizer."""
from itertools import product

import sympy as sp

from verify_finite_torus_slices import conormal, double_middle
from verify_simple_weight_torus import B,C,monomial_chart,rotate,solve_chart

K = range(4)


def quotient(w,weights,top,perm,scales,h,check_sign=False):
    words = [q for q in product(K,repeat=4) if sum(weights[v] for v in q)==top+h]
    conditions = conormal(w,(0,1,2,3))+conormal(w,(1,2,3,0))
    derivative = sp.Matrix([[r.get(q,0) for q in words] for r in conditions if next(iter(r)) in words])
    orbit_rows = [[w.get(q[:i]+(b,)+q[i+1:],0) if q[i]==a else 0 for q in words]
                  for i,a,b in product(K,K,K) if weights[a]-weights[b]==h]
    orbit = sp.Matrix(orbit_rows).T.columnspace() if orbit_rows else []
    columns = sp.Matrix.hstack(*(orbit+derivative.nullspace()))
    basis = columns[:,list(columns.rref()[1])]
    quotient_basis = basis[:,len(orbit):]
    action = sp.Matrix(len(words),quotient_basis.cols,lambda i,j:
                       sp.prod(scales[k][v] for k,v in enumerate(words[i]))*quotient_basis[words.index(perm(words[i])),j])
    shift = basis.gauss_jordan_solve(action)[0][len(orbit):,:]
    if not check_sign:
        return shift,None
    signed = sp.Matrix(len(words),quotient_basis.cols,lambda i,j:(-1)**words[i].count(2)*quotient_basis[i,j])
    sign = basis.gauss_jordan_solve(signed)[0][len(orbit):,:]
    return shift,sign


def verify():
    chart,_,_ = monomial_chart()
    main,line = solve_chart(chart)
    half = {q:v.subs(B,2) for q,v in line.items()}
    cross = {q:v.subs({B:sp.Rational(1,2),C:1}) for q,v in main.items()}
    middle = {q:sp.sympify(v) for q,v in double_middle(sp.Integer(2)).items()}
    half_shift = lambda q:rotate(q,2)
    shift = lambda q:rotate(q,1)
    middle_shift = lambda q:rotate(tuple((0,2,1,3)[v] for v in q),1)
    unit = [1,1,1,1]
    k = sp.sqrt(2)
    half_scales = [[sp.Rational(1,2),1,1,1],unit,[2,1,1,1],unit]
    cross_scales = [unit,unit,[sp.Rational(1,2),1,1,1],[2,1,1,1]]
    middle_scales = [[sp.Rational(1,2),-k/4,sp.Rational(1,4),-sp.Rational(1,8)],
                     [4,-2*k,-2,1],[2,-k,-2,1],[-2,-k/2,1,1]]
    for w,g,d in ((half,half_shift,half_scales),(cross,shift,cross_scales),(middle,middle_shift,middle_scales)):
        for q in product(K,repeat=4):
            assert sp.simplify(sp.prod(d[i][v] for i,v in enumerate(q))*w.get(g(q),0)-w.get(q,0))==0
    for h in range(1,6):
        a,_ = quotient(half,tuple(K),6,half_shift,half_scales,h)
        assert a == sp.eye(a.rows)
        b,_ = quotient(cross,tuple(K),6,shift,cross_scales,h)
        assert b*b == sp.eye(b.rows)
        assert (b-sp.eye(b.rows)).rank()==int(h in (3,4))
    for h in (1,2,3):
        f,delta = quotient(middle,(0,1,1,2),4,middle_shift,middle_scales,h,True)
        assert (f*f-(k/4)**h*sp.eye(f.rows)).applyfunc(sp.simplify).is_zero_matrix
        assert f*delta == (-1)**h*delta*f
        if h==2:
            assert sp.simplify(f.trace()-k/4)==0 and sp.simplify((f*delta).trace()+k/4)==0


if __name__ == '__main__':
    verify()
    print('Half shift fixes the curve slice. The crossing and double-middle shifts act by the specified involutions.')
