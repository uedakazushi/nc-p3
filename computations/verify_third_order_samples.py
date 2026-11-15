#!/usr/bin/env python3
"""Three exact third-order test directions; NOT universal integrability."""
from pbw_complex import *
from verify_moduli_jets import quotient_basis, second_order_obstruction
N,M=matrices(clifford0())
O=quotient_basis(N,M,-1)
L=[s.Matrix.hstack(*(M[0]+z*M[1]).T.nullspace()).T for z in [1,-1]]
rules=clifford0();zero=s.zeros(60,1)
for u in [[1,1,-1,0,0,1],[1,2,3,4,5,-14],[0,1,2,1,0,0]]:
    # Q = u0*u5+u1*u2-u2*u3+u3*u4
    Q=u[0]*u[5]+u[1]*u[2]-u[2]*u[3]+u[3]*u[4]
    assert Q == 0
    alpha=O*s.Matrix(u);ob=second_order_obstruction(rules,zero,alpha)[0]
    sol,pars=(M[0]+M[1]).gauss_jordan_solve(-ob)
    beta=sol.subs({p:0 for p in pars})
    ars={}
    for order in [1,2]:
        v=alpha if order==1 else beta
        ars[order]={tuple(reversed(ab)):{m:v[10*r+j] for j,m in enumerate(MON2) if v[10*r+j]} for r,ab in enumerate(PAIRS)}
    @lru_cache(None)
    def nf(w,start,order):
        for j in range(len(w)-1):
            pair=w[j:j+2]
            if pair[0]>pair[1]:
                d={}
                for m,c in rules[pair].items():
                    for n,v in nf(w[:j]+m+w[j+2:],start,order).items():add(d,n,c*v)
                for deg in range(1,min(order,2)+1):
                    sig=(-1)**(start+j) if deg==1 else 1
                    for m,c in ars[deg][pair].items():
                        for n,v in nf(w[:j]+m+w[j+2:],start,order-deg).items():add(d,n,sig*c*v)
                return d
        return {w:s.Integer(1)} if order==0 else {}
    def forced(w,start,j,order):
        d={}
        for m,c in rules[w[j:j+2]].items():
            for n,v in nf(w[:j]+m+w[j+2:],start,order).items():add(d,n,c*v)
        for deg in range(1,min(order,2)+1):
            sig=(-1)**(start+j) if deg==1 else 1
            for m,c in ars[deg][w[j:j+2]].items():
                for n,v in nf(w[:j]+m+w[j+2:],start,order-deg).items():add(d,n,sig*c*v)
        return d
    vals=[]
    for start in [0,1]:
        out=s.zeros(80,1)
        for r,t in enumerate(TRIPLES):
            w=tuple(reversed(t));d=forced(w,start,0,3)
            for n,v in forced(w,start,1,3).items():add(d,n,-v)
            for n,v in d.items():out[20*r+MON3.index(n)]=v
        vals.append(out)
    proj=L[1]*((vals[0]-vals[1])/2)
    assert proj == s.zeros(L[1].rows,1)
    assert vals[0]+vals[1] == s.zeros(80,1)
    print('PASS: third-order lift for sample', u, flush=True)
