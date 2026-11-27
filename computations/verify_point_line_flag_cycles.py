"""Polynomial Grassmannian incidence certificates for the two rank-ten sheets."""
from itertools import combinations,product,permutations
from functools import reduce
import sympy as sp
from flint import fmpz_mpoly_ctx
from verify_full_thin_exchange_cover import tensor
K=range(4)
def prim(v):
    v=list(v);w=[a for a in v if a]
    if not w: return v
    f=reduce(lambda a,b:a.gcd(b) if not a.is_one() else a,w)
    return [a//f for a in v]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def trans(D):return list(map(list,zip(*D)))
def wedge(a,b):return [a[i]*b[j]-a[j]*b[i] for i,j in combinations(K,2)]
def three(a,b,c):
    return [a[i]*(b[j]*c[k]-b[k]*c[j])-a[j]*(b[i]*c[k]-b[k]*c[i])+a[k]*(b[i]*c[j]-b[j]*c[i]) for i,j,k in combinations(K,3)]
def perp3(a,b,c):
    t=three(a,b,c);return prim([t[3],-t[2],t[1],-t[0]])
def plane(rows):
    a=next(r for r in rows if any(r))
    b=next(r for r in rows if any(wedge(a,r)))
    assert all(not any(three(a,b,r)) for r in rows)
    return [prim(a),prim(b)]
def equal(U,V):return all(not any(three(U[0],U[1],r)) for r in V)
def perp(U):
    a,b=U;w=wedge(a,b)
    i,j=next(ij for ij,v in zip(combinations(K,2),w) if v)
    M=a[i]*b[j]-a[j]*b[i]
    V=[]
    for k in K:
        if k in (i,j):continue
        q=[0*M for _ in K];q[k]=M
        q[i]=a[j]*b[k]-a[k]*b[j];q[j]=a[k]*b[i]-a[i]*b[k]
        V.append(prim(q))
    return V
def project(ctx,vs,base):
    p,vals=base
    ev=[int(v(*vals))%p for v in vs]
    c=next(z for z in ev if z)
    return tuple(z*pow(c,-1,p)%p for z in ev)
def two(T,start,x,y):
    p=start%4;q=(p+1)%4;k=(p+2)%4;l=(p+3)%4
    v=[[0*T[next(iter(T))] for _ in K] for _ in K]
    for w,c in T.items():
        v[w[k]][w[l]]+=c*x[w[p]]*y[w[q]]
    return v
def step(T,p,U,prev,pt):
    assert all(not dot(v,pt) for v in U)
    q=next(v for v in perp(U) if any(wedge(v,pt)))
    return plane(trans(two(T,p-1,prev,q)))
def nextpoint(T,p,x,forward=True):
    p=p%4;q=(p+(1 if forward else -1))%4
    r=(p+(2 if forward else -2))%4;l=(p+(3 if forward else -3))%4
    rows=[[0*T[next(iter(T))] for _ in K] for _ in range(16)]
    for w,c in T.items():rows[4*w[r]+w[l]][w[q]]+=c*x[w[p]]
    a=next(v for v in rows if any(v));b=next(v for v in rows if any(wedge(a,v)));c=next(v for v in rows if any(three(a,b,v)))
    ans=perp3(a,b,c)
    assert all(not dot(z,ans) for z in rows)
    return ans
def pol(ctx,expr,syms):
    return ctx.from_dict({k:int(v) for k,v in sp.Poly(expr,*syms).terms()})

def second_sheet():
    syms=a,b0,c0,d=sp.symbols('a b c d')
    b=b0*b0;c=c0*c0;A=sp.Matrix([[a,b],[c,d]]);B=sp.Matrix([[-d,b],[c,-a]])
    T0=tensor(A,B,-sp.S.One,sp.S.One);ctx=fmpz_mpoly_ctx.get(('a','b','c','d'))
    T={w:pol(ctx,sp.cancel(z*(a*d-b*c)),syms) for w,z in T0.items()}
    a,B0,C0,d=ctx.gens();b=B0**2;c=C0**2;one=ctx.constant(1);zero=ctx.constant(0)
    # Interpolate the choice of sign at a single prime; all subsequent checks are polynomial identities.
    base=109,(2,int(sp.sqrt_mod(3,109)),int(sp.sqrt_mod(5,109)),7)
    targets={1:(0,0,1,1),4:(0,1,47,26),9:(1,0,59,10),10:(1,0,50,99),15:(0,1,62,83)}
    pts={}
    for i,target in targets.items():
        if i==1:vs=[zero,zero,one,one]
        else:
            basis=[zero,2*C0,d+c,d-c] if target[0]==0 else [2*B0,zero,a+b,a-b]
            vs=next([a*sign for a,sign in zip(basis,signs)] for signs in product((1,-1),repeat=4)
                    if project(ctx,[a*sign for a,sign in zip(basis,signs)],base)==target)
        pt=[None,vs,None,None]
        for p in (1,2,3):pt[(p+1)%4]=nextpoint(T,p,pt[p])
        assert not any(wedge(vs,nextpoint(T,0,pt[0])))
        pts[i]=pt
    def pt(i,k):return pts[i][k%4]
    def rel(k,i,j):return two(T,k-2,pt(i,k-2),pt(j,k-1))
    r=[rel(0,4,1),rel(1,1,4),rel(2,4,1),rel(3,1,4)]
    Us=[plane(trans(y)) for y in r]
    # At odd positions the specified contractions give the other odd coordinate plane as left support;
    # the displayed G nodes there are instead obtained as right support on the preceding even node.
    Us[1]=plane(r[0]);Us[3]=plane(r[2])
    for p in (1,3):
        for i,j in permutations(targets,2):
            tt=rel(p,i,j)
            if all(not any(three(Us[p][0],Us[p][1],v)) for v in trans(tt)) and any(v for row in tt for v in row):
                if equal(Us[(p+1)%4],plane(tt)):r[p]=tt;break
        else:raise AssertionError(('relation',p))
    for p in K:assert equal(Us[p],plane(trans(r[p]))) and equal(Us[(p+1)%4],plane(r[p]))

    # Only choose a point for each transition; the proper incidence uses the two-dimensional quotient.
    for p,walk,repeat in [(0,[1,9,10,15,4],1),(1,[15,4,10,9],0),(2,[1,9,10,15,4],1),(3,[15,4,9,10],0)]:
        U=Us[p];past=[]
        for j,i in enumerate(walk):
            past.append(U)
            U=step(T,p+j,U,pt(i,p+j-1),pt(i,p+j))
        assert equal(U,past[repeat])



def conic_sheet():
    syms=d0,k0,l=sp.symbols('d k l')
    b=sp.S.One;d=d0;k=k0;e=sp.factor((k*d-1)/(d-k))
    a=sp.factor(((e-1)-l*l*(e+1))/((e-1)+l*l*(e+1)))
    h=sp.factor((k*a-1)/(a-k));D=sp.factor(a*d-1)
    mu=[sp.factor(2*(a+e)/((a+1)*(e+1))),sp.factor(2*(d+h)/((d+1)*(h+1)))]
    L=[l,-(k-b)/((k+b)*l)]
    ctx=fmpz_mpoly_ctx.get(('d','k','l'))
    def clear(v):
        vs=[sp.factor(t) for t in v];den=sp.lcm_list([sp.fraction(t)[1] for t in vs])
        return prim([pol(ctx,sp.cancel(t*den),syms) for t in vs])
    def rank(rows):
        rows=[prim(r) for r in rows];r=0
        for j in range(len(rows[0])):
            k=next((i for i in range(r,len(rows)) if rows[i][j]),None)
            if k is None:continue
            rows[r],rows[k]=rows[k],rows[r]
            for i in range(r+1,len(rows)):
                if rows[i][j]:
                    rows[i]=prim([y*rows[r][j]-z*rows[i][j] for y,z in zip(rows[i],rows[r])])
            r+=1
        return r
    aa,bb,dd,ee,hh=var=sp.symbols('a b D e h')
    raw=tensor(sp.Matrix([[aa,bb],[bb,dd]]),sp.Matrix([[ee,bb],[bb,hh]]),-sp.S.One,sp.S.One)
    for typ in range(2):
        sc=mu[typ];sub=dict(zip(var,[sp.factor(t*sc) for t in (a,b,d,e,h)]))
        T=dict(zip(raw,clear([v.subs(sub) for v in raw.values()])))
        dc=sp.factor(sc*D);dg=sp.factor(sc*(e*h-1))
        m,n,f,g=[sp.factor(z/dc) for z in (d-1,d+1,a-1,a+1)]
        M,N,F,G=[sp.factor(z/dg) for z in (h-1,h+1,e-1,e+1)]
        sec={}
        for t in (0,1):
            L0=L[t]
            for p in K:
                if p%2==0:
                    rows=[[1-t,t,0,0],[0,0,1,L0]]
                else:
                    z=L0 if p==1 else -L0
                    if t==0:rows=[[-(m*g+n*f)*z,M*f-N*g*z*z,0,0],[0,0,g*z,-f]]
                    else:rows=[[G*n*z*z-F*m,-(f*n+g*m)*z,0,0],[0,0,n*z,m]]
                sec[p,t]=perp([clear(v) for v in rows])
        pts=[]
        for sign in (1,-1):
            pt=[clear([1-typ,typ,sign,sign*L[typ]])]
            for p in range(3):pt.append(nextpoint(T,p,pt[-1]))
            assert not any(wedge(pt[0],nextpoint(T,3,pt[3])))
            for p in K:assert all(not dot(row,pt[p]) for row in sec[p,typ])
            pts.append(pt)
        for p in K:
            # Each secant maps by two successive point steps to the other secant two positions later.
            W=step(T,p,sec[p,typ],pts[0][(p-1)%4],pts[0][p])
            WW=step(T,p+1,W,pts[1][p],pts[1][(p+1)%4])
            assert equal(WW,sec[(p+2)%4,1-typ])
        if typ==0:
            for p in K:
                positions=[(p+i)%4 for i in K]
                def rv(i,j,k,n):
                    q=[None]*4
                    for pos,c in zip(positions,(i,j,k,n)):q[pos]=c
                    return T.get(tuple(q),ctx.constant(0))
                rows=[[rv(i,j,k,n) for i,j in product(K,repeat=2)] for k,n in product(K,repeat=2)]
                extra=[[i*j for i,j in product(u,v)] for u,v in product(sec[p,p%2],sec[(p+1)%4,1-p%2])]
                assert rank(rows)==6 and rank(rows+extra)==9



if __name__=='__main__':
    second_sheet()
    conic_sheet()
    print('PASS: generic rank-two relation cycles and all recurrent point steps on both sheets')
