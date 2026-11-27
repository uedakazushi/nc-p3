"""Unrestricted normal weights and four-slot evaluations in three Artin point algebras."""
from itertools import product

import sympy as sp

from verify_double_middle_torus import cubics, forms, recurrence
from verify_simple_weight_torus import B, C, monomial_chart, rotate, solve_chart

K = range(4)
S,T = sp.symbols('S T')


def double_middle(a):
    low,high = forms('diagonal',a,1/a)
    mixed = recurrence(low,high,sp.eye(2))
    w = {(0,3,0,3):1,(3,0,3,0):1,(1,1,1,1):1,(2,2,2,2):1}
    for i in K:
        for j,k in product((1,2),repeat=2):
            for word,mat in (((0,3,j,k),low[i]),((3,0,j,k),high[i]),((0,j,3,k),mixed[i])):
                w[rotate(word,(-i)%4)] = mat[j-1,k-1]
    return w


def conormal(w, order):
    def word(i,j):
        q = (i//4,i%4,j//4,j%4)
        return tuple(q[order.index(k)] for k in K)
    cut = sp.Matrix(16,16,lambda i,j:w.get(word(i,j),0))
    return [{word(i,j):l[i]*r[j] for i in range(16) if l[i] for j in range(16) if r[j]}
            for l,r in product(cut.T.nullspace(),cut.nullspace())]


def normal_table(w,weights,top,sign=False):
    rows = conormal(w,(0,1,2,3))+conormal(w,(1,2,3,0))
    table = []
    for total in range(4*max(weights)+1):
        for parity in range(2 if sign else 1):
            words = [q for q in product(K,repeat=4) if sum(weights[v] for v in q)==total
                     and (not sign or q.count(2)%2==parity)]
            eq = [[r.get(q,0) for q in words] for r in rows if next(iter(r)) in words]
            orbit = [[w.get(q[:i]+(b,)+q[i+1:],0) if q[i]==a else 0 for q in words]
                     for i,a,b in product(K,K,K) if weights[a]-weights[b]==total-top
                     and (not sign or int((a==2) != (b==2))==parity)]
            table.append(len(words)-sp.Matrix(eq).rank()-sp.Matrix(orbit).rank())
    return table


def simple_local(w,rev,start):
    t,z = S,T
    slots = []
    for step in K:
        i = (step+start)%4
        def val(*q):
            if rev:
                q = tuple(3-v for v in q)
            return w.get(rotate(q,(-i)%4),0)
        v = -val(1,0,2,3)*t/val(0,1,2,3)
        zz = -(val(1,1,1,3)*t*v+val(2,0,1,3)*z)/val(0,2,1,3)
        x3 = -(val(1,2,3,0)*t*zz+val(2,1,3,0)*z*v)/val(3,0,3,0)
        y3 = -(val(1,2,0,3)*t*zz+val(2,1,0,3)*z*v)/val(0,3,0,3)
        x,y = (sp.Integer(1),t,z,sp.expand(x3)),(1,v,zz,y3)
        if step==0:
            eq = [sum(val(a,n-a,k,l)*x[a]*y[n-a] for a in K if 0<=n-a<4)
                  for n,(k,l) in [(4,(0,2)),(5,(0,1))]]
            gb = sp.groebner(eq,S,T,order='grlex',domain=sp.QQ)
        slots.append(x[::-1] if rev else x)
        t,z = sp.expand(v),sp.expand(zz)
    return slots,gb


def middle_local(rev,start):
    low,high = forms('diagonal',sp.Integer(2),sp.Rational(1,2))
    mixed = recurrence(low,high,sp.eye(2))[:4]
    if rev:
        low,high,mixed = high,low,[mixed[(i+2)%4].T for i in K]
    u = sp.Matrix([S,T])
    slots = []
    poly = cubics(low,high,mixed,start)
    gb = sp.groebner([e.subs(dict(zip(sp.symbols('X Y'),u))) for e in poly],S,T,order='grlex',domain=sp.QQ)
    for step in K:
        i = (step+start)%4
        v = -mixed[i].inv().T*low[(i+1)%4]*u
        z = -(u.T*high[(i+2)%4]*v)[0]
        slots.append((z,*u,1) if rev else (1,*u,z))
        u = v
    return slots,gb


def evaluation_rank(local,middle=False):
    rows = []
    for rev in (False,True):
        slots,gb = local(rev,3)
        values = [sp.Poly(gb.reduce(sp.expand(sp.prod(slots[i][q[i]] for i in K)))[1],S,T).as_dict()
                  for q in product(K,repeat=4)]
        for degree in sorted(set().union(*(p.keys() for p in values))):
            rows.append([v.get(degree,0) for v in values])
    if middle:
        for j in (0,1):
            rows.append([int(all(q[i]==1+(i+j)%2 for i in K)) for q in product(K,repeat=4)])
    return sp.Matrix(rows).rank()


def verify():
    chart,_,_ = monomial_chart()
    main,line = solve_chart(chart)
    main = {w:v.subs({B:2,C:3}) for w,v in main.items()}
    line = {w:v.subs(B,2) for w,v in line.items()}
    middle = double_middle(sp.Integer(2))
    for w,center,ev in ((main,2,20),(line,1,19)):
        assert normal_table(w,(0,1,2,3),6)==[0,1,2,2,1,1,center,1,1,2,2,1,0]
        assert evaluation_rank(lambda r,i:simple_local(w,r,i)) == ev
    assert normal_table(middle,(0,1,1,2),4)==[0,2,3,2,1,2,3,2,0]
    assert normal_table(middle,(0,1,1,2),4,True)[10:16]==[1,1,2,1,1,1]
    assert evaluation_rank(middle_local,True)==19


if __name__ == '__main__':
    verify()
    print('Unrestricted smooth orbit-normal slices: 16,15,15; four-slot evaluation ranks: 20,19,19.')
