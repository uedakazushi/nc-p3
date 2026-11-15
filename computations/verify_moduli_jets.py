#!/usr/bin/env python3
"""Reproduce the nc-p3 moduli tangent and quadratic-obstruction certificates.

Run from the repository root:
    python3 computations/verify_moduli_jets.py --output research/moduli_certificate.json
This is an exact finite computation, not formal verification of the entire
mathematical argument, and not a claim of higher-order unobstructedness.
"""
from pbw_complex import *
from pathlib import Path
import argparse
import json
import platform


def quotient_basis(N, M, z):
    D1=N[0]+z*N[1];D2=M[0]+z*M[1]
    K=s.Matrix.hstack(*D2.nullspace())
    G=s.Matrix.hstack(*D1.columnspace())
    piv=(G.row_join(K)).rref()[1]
    H=K[:,[j-G.cols for j in piv if j>=G.cols]]
    assert D2*H == s.zeros(80,H.cols)
    assert rank(G.row_join(H)) == G.cols+H.cols
    return H


def second_order_obstruction(rules, even, odd):
    """Order-t^2 overlap residual; output its even and odd Fourier modes."""
    ars=[]
    for sig in [1,-1]:
        alpha=even+sig*odd
        ars.append({tuple(reversed(ab)):{m:alpha[10*r+j] for j,m in enumerate(MON2) if alpha[10*r+j]} for r,ab in enumerate(PAIRS)})
    @lru_cache(None)
    def nf(w,start,order):
        for j in range(len(w)-1):
            pair=w[j:j+2]
            if pair[0]>pair[1]:
                d={}
                for m,c in rules[pair].items():
                    for n,v in nf(w[:j]+m+w[j+2:],start,order).items():add(d,n,c*v)
                if order:
                    for m,c in ars[(start+j)%2][pair].items():
                        for n,v in nf(w[:j]+m+w[j+2:],start,order-1).items():add(d,n,c*v)
                return d
        return {w:s.Integer(1)} if order==0 else {}
    def forced(w,start,j):
        d={}
        for m,c in rules[w[j:j+2]].items():
            for n,v in nf(w[:j]+m+w[j+2:],start,2).items():add(d,n,c*v)
        for m,c in ars[(start+j)%2][w[j:j+2]].items():
            for n,v in nf(w[:j]+m+w[j+2:],start,1).items():add(d,n,c*v)
        return d
    out=[]
    for start in [0,1]:
        vout=s.zeros(80,1)
        for r,t in enumerate(TRIPLES):
            w=tuple(reversed(t));d=forced(w,start,0)
            for n,v in forced(w,start,1).items():add(d,n,-v)
            for n,v in d.items():vout[20*r+MON3.index(n)]=v
        out.append(vout)
    return ((out[0]+out[1])/2,(out[0]-out[1])/2)


def verify_axis(pair, mon):
    t=s.symbols('t');rs=[clifford0(),clifford0()]
    for i in range(2):
        rs[i][pair]=dict(rs[i][pair]);rs[i][pair][mon]=rs[i][pair].get(mon,0)+(-1)**i*t
    @lru_cache(None)
    def nf(w,start):
        for j in range(len(w)-1):
            if w[j]>w[j+1]:
                d={}
                for m,c in rs[(start+j)%2][w[j:j+2]].items():
                    for n,v in nf(w[:j]+m+w[j+2:],start).items():add(d,n,s.expand(c*v))
                return d
        return {w:s.Integer(1)}
    for start in [0,1]:
        for tri in TRIPLES:
            w=tuple(reversed(tri));outputs=[]
            for j in [0,1]:
                d={}
                for m,c in rs[(start+j)%2][w[j:j+2]].items():
                    for n,v in nf(w[:j]+m+w[j+2:],start).items():add(d,n,s.expand(c*v))
                outputs.append(d)
            assert outputs[0] == outputs[1]
    return True


def sparse_basis(H):
    return [[{'pair':list(PAIRS[r//10]),'monomial':list(MON2[r%10]),'coefficient':str(H[r,j])}
             for r in range(60) if H[r,j]] for j in range(H.cols)]


def verify():
    data={'date':'2026-11-15','version':8,'arithmetic':'exact QQ and QQ(i)',
          'python':platform.python_version(),'sympy':s.__version__}
    q=[s.Rational(2),s.Rational(3),s.Rational(1,6),s.Rational(5),s.Rational(2,5),s.Rational(15)]
    params=s.symbols('q0:6')
    qq={(a,b):q[j] for j,(a,b) in enumerate(PAIRS)}
    for a in range(4):
        assert s.prod(qq[(a,b)] if a<b else 1/qq[(b,a)] for b in range(4) if b!=a)==1
    for a,b,c in TRIPLES:assert qq[a,b]*qq[b,c]/qq[a,c] != 1
    ranks={}
    for name,rules in [('quantum',quantum(q)),('clifford',clifford0())]:
        N,M=matrices(rules);rows=[]
        for z in [s.Integer(1),s.Integer(-1),s.I,-s.I]:
            r1=rank(N[0]+z*N[1]);r2=rank(M[0]+z*M[1]);h=60-r1-r2
            expected = ([12,42,6] if z==1 else [15,45,0]) if name=='quantum' else ([15,36,9] if z==1 else [15,39,6] if z==-1 else [15,45,0])
            assert [r1,r2,h] == expected
            rows.append({'mode':str(z),'rank_d1':r1,'rank_d2':r2,'h2':h})
        ranks[name]=rows
        if name=='quantum':
            KS=s.zeros(60,6)
            for j,pair in enumerate(PAIRS):KS[10*j+MON2.index(pair),j]=1
            assert (M[0]+M[1])*KS == s.zeros(80,6)
            assert rank((N[0]+N[1]).row_join(KS))==18
            data['quantum']={'parameters':[str(x) for x in q], 'ks_quotient_rank':6,
                             'graded_parameter_dimension':6,'indexed_z_parameter_dimension':3}
    E=quotient_basis(N,M,1);O=quotient_basis(N,M,-1)
    L=[s.Matrix.hstack(*(M[0]+z*M[1]).T.nullspace()).T for z in [1,-1]]
    assert E.cols==9 and O.cols==6
    zero=s.zeros(60,1)
    basis=[(E[:,j],zero) for j in range(9)]+[(zero,O[:,j]) for j in range(6)]
    var=s.symbols('e0:9')+s.symbols('o0:6')
    diag=[];pol=s.zeros(85,1)
    for j,(e,o) in enumerate(basis):
        ob=second_order_obstruction(clifford0(),e,o)
        v=(L[0]*ob[0]).col_join(L[1]*ob[1]);diag.append(v);pol+=var[j]**2*v
    for a,b in combinations(range(15),2):
        ob=second_order_obstruction(clifford0(),basis[a][0]+basis[b][0],basis[a][1]+basis[b][1])
        v=(L[0]*ob[0]).col_join(L[1]*ob[1])-diag[a]-diag[b]
        pol+=var[a]*var[b]*v
    o=var[9:]
    Q=o[0]*o[5]+o[1]*o[2]-o[2]*o[3]+o[3]*o[4]
    rows=[s.factor(v) for v in pol if s.expand(v)!=0]
    assert rows and all(s.expand(v-Q)==0 or s.expand(v+Q)==0 for v in rows)
    assert s.hessian(Q,o).det()!=0
    assert verify_axis((1,0),(2,3)) and verify_axis((3,2),(0,1))
    # The two simple axes are precisely the first and last odd representatives.
    assert O[10*PAIRS.index((0,1))+MON2.index((2,3)),0]==1
    assert sum(bool(x) for x in O[:,0])==1
    assert O[10*PAIRS.index((2,3))+MON2.index((0,1)),5]==1
    assert sum(bool(x) for x in O[:,5])==1
    data.update({'result':'PASS','ranks':ranks,'even_basis':sparse_basis(E),
                 'odd_basis':sparse_basis(O),'quadratic_obstruction':str(Q),
                 'quadratic_form_rank':rank(s.hessian(Q,o)),
                 'full_15_variable_quadratic_obstruction_span_dimension':1,
                 'axis_a_all_polynomial_overlaps':True,'axis_b_all_polynomial_overlaps':True,
                 'limits':['No complete higher-order Kuranishi equation computed.',
                           'No claim that the Clifford local component has dimension 14.',
                           'All choices of a one-period structure have not been classified.',
                           'This code checks finite identities, not the full external as4 periodicity theorem.']})
    family_ranks=[]
    for lam in [1,2]:
        for period in [2,4]:
            rules=[clifford0() for _ in range(period)]
            for i in range(period):rules[i][(3,2)][(0,1)]=(-1)**i*lam
            GN,GM=periodic_matrices(rules,period)
            rn,rm=rank(GN),rank(GM)
            assert [rn,rm,60*period-rn-rm] == ([30,76,14] if period==2 else [60,166,14])
            family_ranks.append({'lambda':lam,'period':period,'rank_d1':rn,'rank_d2':rm,'h2':14})
    data['period_two_family_tangent_checks']=family_ranks
    return data

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=verify()
    text=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['result','ranks','quadratic_obstruction','quadratic_form_rank','limits']},indent=2))
