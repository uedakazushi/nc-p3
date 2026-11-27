"""The three-parameter completion on the remaining singular-corner thin pattern."""
from verify_parallel_invertible_rectangles import K, WORDS, bit, candidates, matrix, product, s, smith
from verify_klein_three_lines import thin_system
from verify_transverse_rank_twelve_isomorphism import alternating_gradient
from flint import nmod_mat


def third_rectangle(A,D,b):
    A,D,b = map(s.sympify,(A,D,b))
    L = A*D-1
    M = (b-1)*L
    C = 1-(b-1)*A*D/b
    u = -1/L
    values = (A,b,1,1,1,1,C,D,u,u,(A*D-b)/(b*L),D/L,A/L,-b/M,u,u)
    mask = candidates(1)[86]
    keys = [k for k in WORDS if (k[0]^k[1]<2 and k[1]^k[2]<2) or mask&bit(k)]
    holes = [k for k in keys if k[0]^k[1]<2 and k[1]^k[2]<2]
    w = dict.fromkeys(keys,s.S.One)
    w.update(zip(holes,values))
    return {k:s.factor(v) for k,v in w.items()}


def check():
    A,D,b = s.symbols('A D b',nonzero=True)
    w = third_rectangle(A,D,b)
    thin,ds,_,H = thin_system([(2,3)]*2,list(w))
    assert all(d==1 for d in ds) and H.rows==H.rank()
    gauge = s.Matrix([[int(k[p]==i) for p,i in product(K,repeat=2)] for k in thin])
    d,_,r = smith(gauge)
    assert [n for n in d.diagonal() if n]==[1]*11+[2]
    holes = [k for k in w if k not in thin]
    action = s.Matrix([[int(k[p]==i) for p,i in product(K,repeat=2)] for k in holes])*r
    assert all(n%2==0 for n in action[:,11]) and action[:,12:].rank()==1
    assert all(action[j,12:]==action[0,12:] for j in range(8))
    for p,h in product(range(2),K):
        m = matrix(w,p,h)
        if h>=2:
            assert m.rank()==1
            continue
        assert m.rank()==2
        if p==1:
            assert all(s.factor(m[a:a+2,b0:b0+2].det())!=0 for a,b0 in ((0,0),(0,2),(2,0),(2,2)))
    # The pivot minors forcing both inverses in the first cut are present on all three masks.
    for i in (88,89,90):
        M = candidates(1)[i]
        assert all(M&bit(holes[j]) for j in (1,13))
    pts = ((1,2,2,88),(1,3,3,89),(2,3,5,90))
    for a,d,b0,i in pts:
        W = third_rectangle(a,d,b0)
        assert sum(bit(k) for k,v in W.items() if v)==candidates(1)[i]
    for vals in ((1,2,2),(2,3,5)):
        sample = third_rectangle(*vals)
        def residue(z):
            a,b0 = s.fraction(z)
            return int(a)*pow(int(b0),-1,1009)%1009
        for p,h in product(K,K):
            m = matrix(sample,p,h)
            assert nmod_mat([[residue(v) for v in row] for row in m.tolist()],1009).rank()==1+(h<2)
            assert any(residue(v) for k,v in sample.items() if k[p]==h)
    sample = {k:sample.get(k,0) for k in WORDS}
    assert all(s.Matrix([[sample[i,j,i^h,j^h] for j in K] for i in K]).det()!=0 for h in K)
    lie = [[residue(sample.get(k[:p]+(j,)+k[p+1:],0)) if i==k[p] else 0
            for p,i,j in product(K,repeat=3)] for k in product(K,repeat=4)]
    assert nmod_mat(lie,1009).rank()==61
    G = alternating_gradient(sample)
    cap = sum(sample[k]*G[k] for k in WORDS)/4
    inv = []
    for p in K:
        blocks = [matrix(sample,p,h)*matrix(G,p,h).T/cap for h in K]
        inv.append(s.factor(sum(sum((m*m)[i,i] for m in blocks)**2 for i in K)))
    assert len(set(inv))==4
    print('PASS: normalized thin torus, all four completions, pivot cells and all three support strata')


if __name__=='__main__':
    check()
