"""Four-start regularity and intrinsic frames on the rank-fourteen surface and its diagonal."""
from flint import nmod_mat

from verify_parallel_curve_residues import residue_test
from verify_parallel_horizontal_moduli import window
from verify_parallel_invertible_rectangles import K, WORDS, alternating_gradient, matrix, s


def check(diagonal):
    P = 1009
    u = s.Integer(2)
    def residue(x):
        a,b = s.fraction(x)
        return int(a)*pow(int(b),-1,P)%P
    X = s.Symbol('X')
    if diagonal:
        v = s.Integer(729)
        A,B = u+1/u,v+1/v
        z = w = s.Integer(residue((A-B)/(A+B+4)))
        f = 4*X**4+18*X**3+37*X**2+18*X+4
        r = v
    else:
        v = s.Integer(3)
        z = s.Integer(883)
        w = s.Integer(residue(-s.Rational(10,59)-z))
        f = 59*X**2+10*X-81
        r = z
    assert s.Poly(f).is_irreducible and f.subs(X,r)%P==0 and f.diff(X).subs(X,r)%P!=0
    C = (u+v)*(1+u*v)
    assert residue(C*(z+w-2)-u*(v+1)**2*(z*w-1))==0
    assert residue(C*(z+w+2)+v*(u+1)**2*(z*w-1))==0
    t = (u+v)/(1+u*v)
    W = {k:s.Integer(residue(x)) for k,x in window(t,u*t,v*t,u*(z+v)/(1+u*v),u*(w+v)/(1+u*v)).items()}
    assert len(W)==52 and all(W.values())
    raw = {k:W.get(k,0) for k in WORDS}
    G = alternating_gradient(raw)
    assert sum(raw[k]*G[k] for k in WORDS)%P!=0
    rel = []
    for p in K:
        rules = []
        moment = [0]*4
        for h in K:
            M = matrix(W,p,h)
            reduced,rank = nmod_mat([[int(x)%P for x in row] for row in M.tolist()],P).rref()
            assert rank==(2,2,1,1)[h]
            pivots = [next(j for j,x in enumerate(row) if x) for row in reduced.tolist()[:rank]]
            rules.extend((h,{(i,i^h):int(M[i,j]) for i in K if M[i,j]}) for j in pivots)
            operator = M*matrix(G,p,h).T
            moment = [a+operator.row(i).dot(operator.col(i)) for i,a in enumerate(moment)]
        assert len({m%P for m in moment})==4
        rel.append(rules)
    residue_test(0,P,lambda _,prime:rel)


if __name__=='__main__':
    check(False)
    check(True)
    print('PASS: four regular starts and four distinct eigenlines at both algebraic sample points')
