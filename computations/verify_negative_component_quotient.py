"""Thin gauge, projective parameter transformations, and the finite quotient of the negative chart."""
from itertools import product
from functools import reduce
from operator import xor
import sympy as s
from flint import fmpz_mat
from verify_klein_elliptic_completion import gram_twist

K = range(4)
a,b,c,d,z = params = s.symbols('a b c d z')
D = a*d-b*c
words = [q for q in product(K,repeat=4) if reduce(xor,q)==0]
thin = [q for q in words if q[0]^q[1]>=2 or q[1]^q[2]>=2]


def coefficients(A,z):
    w,_,_ = gram_twist(A,-A,z,s.S.NegativeOne)
    return {q:s.factor(w[4*q[0]+q[1],4*q[2]+q[3]]) for q in words}


def lattice(keys):
    E = fmpz_mat([[int(q[p]==i) for p in K for i in K] for q in keys])
    H,T = E.hnf(transform=True)
    assert H==T*E and abs(int(T.det()))==1
    return fmpz_mat(T.tolist()[E.rank():])


def relabel(q,shift,shear,offsets):
    v = tuple((j^(j>>1) if shear else j)^i for j,i in zip(q,offsets))
    return v[shift:]+v[:shift]


def verify_thin(F):
    ker = lattice(thin)
    def weights(fs):
        out = []
        for v in fs:
            u,e = v.as_coeff_exponent(z)
            assert u in (1,-1)
            out.append([int(u==-1),int(e)])
        return fmpz_mat(out)
    for shift,shear,i in product(K,range(2),K):
        for j,k in product((i&2,(i&2)^1),repeat=2):
            off = i,j,k,i^j^k
            G = {q:F[relabel(q,shift,shear,off)] for q in thin}
            w = s.factor(G[0,0,2,2]*G[1,1,3,3]/(G[0,0,3,3]*G[1,1,2,2]))
            assert w in (z,1/z)
            test = ker*(weights([G[q] for q in thin])-weights([F[q].subs(z,w) for q in thin]))
            eq = [z**abs(int(e))-(-1)**(int(t)%2) for t,e in test.tolist() if e or int(t)%2]
            assert (s.gcd_list(eq).expand() if eq else 0)==(0 if i^k==shift%2 else z-1)


def verify_generators(F):
    ker = lattice(words)
    gens = [(0,0,(0,1,0,1),(a,c,b,d),1/z),
            (0,0,(1,0,1,0),(d,z*b,c/z,a),1/z),
            (0,1,(0,0,0,0),(a,z*b,c/z,d),1/z),
            (0,1,(2,2,2,2),(z*(a+b+c+d),-(a*z+b*z+c+d),-z*(a+b*z+c+d*z),a*z+b*z*z+c+d*z),z),
            (1,0,(0,0,1,1),(c,a,z*d,z*b),z)]
    for shift,shear,off,newA,newZ in gens:
        G = {q:s.factor(v.xreplace(dict(zip(params,(*newA,newZ))))) for q,v in F.items()}
        factors = []
        for q in words:
            r = s.factor(G[q]/F[relabel(q,shift,shear,off)])
            unit,fact = s.factor_list(s.fraction(r)[0])
            unit2,den = s.factor_list(s.fraction(r)[1])
            assert abs(unit/unit2)==1
            ex = dict(fact)
            for v,e in den:
                ex[v] = ex.get(v,0)-e
            assert set(ex)<={z,z-1,D}
            factors.append([int(unit/unit2==-1),*[int(ex.get(v,0)) for v in (z,z-1,D)]])
        result = ker*fmpz_mat(factors)
        assert all(int(row[0])%2==0 and not any(row[1:]) for row in result.tolist())
    # With the thin frame fixed, the only thick freedom is a common multiplier of A.
    thinrank = lattice(thin).nrows()
    assert len(thin)-thinrank==12
    extended = thin+[q for q in words if max(q)<2]
    T = lattice(extended)
    top = {(0,0,0,0):a,(0,0,1,1):b,(1,1,0,0):c,(1,1,1,1):d}
    tag = fmpz_mat([[int(F[q] in (v,-v)) for v in top.values()] for q in extended])
    image = (T*tag).hnf()
    assert image.tolist()[:3]==[[1,0,0,-1],[0,1,0,-1],[0,0,1,-1]]
    # This test also covers zeros: the ratios between the four coordinate multipliers
    # are trivial on the full kernel of the thin action.
    rows = [[int(q[p]==i) for p in K for i in K] for q in thin]
    basis = fmpz_mat(rows).hnf()
    matrices = [fmpz_mat(rows+[[int(q[p]==i) for p in K for i in K]]).hnf() for q in top]
    full = fmpz_mat([[int(q[p]==i) for p in K for i in K] for q in words]).hnf()
    assert all(m.tolist()[:13]==full.tolist()[:13] for m in matrices)
    for q in top:
        diff = [int(q[p]==i)-int(i==0) for p in K for i in K]
        augmented = fmpz_mat(rows+[diff]).hnf()
        assert augmented.tolist()[:12]==basis.tolist()[:12]


def quotient():
    t,Z,u,v = s.symbols('t Z u v')
    t1 = -(1+2*t)/(2+Z*t)
    delta = (Z-4)/(2+Z*t)**2
    tr = {t:t1,u:delta*v,v:delta*u}
    for f in (t+t1,(1+delta)*(u+v),(1+delta)*(u-v)/(t-t1)):
        assert s.factor(f.xreplace(tr)-f)==0


if __name__=='__main__':
    F = coefficients(s.Matrix([[a,b],[c,d]]),z)
    verify_thin(F)
    verify_generators(F)
    quotient()
    print('PASS: exact thin exceptions, all five parameter transformations and the rational quotient')
