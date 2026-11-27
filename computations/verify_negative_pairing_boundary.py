"""Compute both canonical boundary frames, their integral gauge lattices, and all support maps."""
from collections import Counter
from itertools import product
import sympy as s
from flint import fmpz_mat

from verify_negative_component_quotient import D,K,a,b,c,d,z,coefficients,lattice
from verify_negative_component_point_graph import evaluation,x
from verify_klein_elliptic_completion import gram_twist

w,m = s.symbols('w m')
had = s.Matrix([[1,1],[1,-1]])


def frame(tensor,bases):
    out = {}
    for q in product(K,repeat=4):
        val = s.factor(sum(f*s.prod(B[i,j] for B,i,j in zip(bases,k,q))
                           for k,f in tensor.items() if all(B[i,j] for B,i,j in zip(bases,k,q))))
        if val:
            out[q] = val
    return out


def action(h,q):
    return tuple((i^((flag>>(i>>1))&1))^(2*(flag>>2)) for i,flag in zip(q,h))


def exponent_matrix(keys):
    return fmpz_mat([[int(q[p]==i) for p in K for i in K] for q in keys])


def invariants(F,parameters,keys):
    ker = lattice(keys)
    rows = []
    for q in keys:
        nu,de = s.fraction(F[q])
        facts = dict(s.factor_list(nu)[1])
        for f,e in s.factor_list(de)[1]:
            facts[f] = facts.get(f,0)-e
        assert set(facts)<=set(parameters)
        rows.append([int(facts.get(v,0)) for v in parameters])
    H,T = (ker*fmpz_mat(rows)).hnf(transform=True)
    return H.tolist()[:H.rank()],(T*ker).tolist()[:H.rank()]


def moved(F,keys,rows,g):
    return tuple(s.factor(s.prod(F[action(g,q)]**int(i) for q,i in zip(keys,row) if i)) for row in rows)


def single():
    F = frame(coefficients(s.Matrix([[a,b],[c,a]]),z),[s.diag(had,s.eye(2)),s.diag(s.eye(2),had)]*2)
    const = [q for q,f in F.items() if f.free_symbols<={z}]
    models = [next(q for q,f in F.items() if f==2*v) for v in (a,b,c)]
    top = const+models
    H,rows = invariants(F,(a,b,c,z),top)
    assert H==[[1,0,-1,0],[0,1,-1,0],[0,0,0,1]]
    norm = tuple(s.factor(f/v) for f,v in zip(moved(F,top,rows,(0,0,0,0)),(a/c,b/c,z)))
    assert len(F)==58 and len(const)==42
    assert exponent_matrix(const).rank()==11
    for q in models:
        mat = exponent_matrix(const+[q])
        assert mat.rank()==12 and [mat.snf()[j,j] for j in range(12)]==[1]*12
        assert mat.hnf().tolist()[:12]==exponent_matrix(F).hnf().tolist()[:12]
    allmaps = Counter()
    failures = 0
    for g in product(K,repeat=4):
        if any(action(g,q) not in F for q in const):
            continue
        assert all(action(g,q) in const for q in const)
        assert all(action(g,q) in F for q in F)
        fs = tuple(s.factor(v/k) for v,k in zip(moved(F,top,rows,g),norm))
        if s.factor(fs[1]*fs[2]-1)==0:
            # Both b' and c' are carried by the same source coordinate, also at its zeros.
            assert s.factor(F[action(g,models[1])]/F[action(g,models[2])]) in (1,-1)
            failures += 1
        else:
            allmaps[fs] += 1
    assert failures==8
    assert allmaps==Counter({(a/c,b/c,z):2,(a/b,c/b,1/z):2,
                              (a*z/c,b*z*z/c,1/z):2,(a/(b*z),c/(b*z*z),z):2})


def double():
    wh,hw = s.diag(w,1)*had,s.diag(1,w)*had
    bases = [s.diag(had,wh),s.diag(wh,had),s.diag(had,wh),s.diag(hw,had)]
    F = frame(coefficients(s.Matrix([[m*w,1],[w*w,m*w]]),w*w),bases)
    keys = list(F)
    H,rows = invariants(F,(w,w+1,w-1,m+1,m-1),keys)
    assert H==[[1,0,-2,0,0],[0,2,-2,0,0],[0,0,0,1,-1]]
    X,Y = (w+1)**2/(w-1)**2,(m+1)/(m-1)
    v0 = moved(F,keys,rows,(0,0,0,0))
    norm = [s.factor(f/g) for f,g in zip(v0,(X-1,X,Y))]
    counts = Counter()
    for g in product(range(8),repeat=4):
        if all(action(g,q) in F for q in keys):
            X0,Y0 = [s.factor(f/n) for f,n in zip(moved(F,keys,rows,g)[1:],norm[1:])]
            counts[X0,Y0] += 1
    assert counts==Counter({(s.factor(u),s.factor(v)):2 for u,v in
                            [(X,Y),(1/X,1/Y),(X,-1/(X*Y)),(1/X,-X*Y)]})
    E = exponent_matrix(keys)
    assert len(keys)==28 and E.rank()==11 and [E.snf()[i,i] for i in range(11)]==[1]*11


def local_minors():
    A = a*a-b*c
    E = evaluation().subs(d,a)
    L = s.symbols('l:4')
    S,T = 2*a+b+c,2*a*z+b*z*z+c
    Q = a*(L[0]-L[1])-L[2]+z*L[3]
    sub = {L[0]:L[1]+(S*L[2]-T*L[3])/A}
    for rows in ('0123','0245','0345','1245','1345'):
        f = s.Poly(E.extract(tuple(map(int,rows)),K).det(method='domain-ge'),*x)
        assert all(all(v%2==0 for v in ex) for ex in f.monoms())
        h = sum(v*s.prod(y**(i//2) for y,i in zip(L,ex)) for ex,v in f.terms())
        g = s.Poly(s.factor(h.subs(sub)/Q.subs(sub)),*L)
        assert all(sum(ex)==1 and ex[:2]==(0,0) for ex in g.monoms())
    W,_,_ = gram_twist(s.Matrix([[a,b],[c,d]]),-s.Matrix([[a,b],[c,d]]),z,s.S.NegativeOne)
    p,r = a+d,z*b+c
    expected = [-4*p*(2*z*p+(z+1)*r)/D,4*r*((z+1)*p+2*r)/D,-4*(z+1)**2,-4*(z+1)**2]
    for h in K:
        block = s.Matrix(4,4,lambda i,j:W[4*i+j,4*(i^h)+(j^h)])
        assert s.factor(block.det()-expected[h])==0


if __name__=='__main__':
    single()
    double()
    local_minors()
    print('PASS: all boundary support permutations, parameter actions, double-conic torsion and opposite determinants')
