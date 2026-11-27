"""The projective vertical boundary, including its allowed shears and linear quotient."""
from itertools import product
from collections import defaultdict
import sympy as s
from flint import fmpz_mat

from verify_negative_component_quotient import K,a,b,c,d,coefficients,lattice
from verify_negative_pairing_boundary import action,exponent_matrix,frame,had
from verify_parallel_invertible_rectangles import smith
from verify_third_rectangle_orbits import support_maps


def generators():
    F = {q:v for q,v in coefficients(s.Matrix([[a,b],[c,d]]),s.S.Zero).items() if v}
    D = a*d-b*c
    for shift,new in ((0,(d,b,c,a)),(1,(-a-c,b+a+c+d,c,-c-d))):
        G = coefficients(s.Matrix(2,2,new),s.S.Zero)
        rotated = {q[shift:]+q[:shift]:f for q,f in F.items()}
        passed = False
        for ps in support_maps(rotated,F):
            src = {tuple(p[i] for p,i in zip(ps,q)):f for q,f in rotated.items()}
            factors = [s.factor(G[q]/src[q]) for q in F]
            units = {x*D**j:(int(x<0),j) for x in (-1,1) for j in (-1,0,1)}
            if all(f in units for f in factors):
                res = lattice(list(F))*fmpz_mat([list(units[f]) for f in factors])
                if all(int(r[0])%2==0 and r[1]==0 for r in res.tolist()):
                    passed = True
        assert passed


def flags(equal):
    F = {q:v for q,v in coefficients(s.Matrix([[a,b],[0,a if equal else d]]),s.S.Zero).items() if v}
    if equal:
        F = frame(F,[s.diag(had,s.eye(2)),s.diag(s.eye(2),had)]*2)
    vertex,lift = (3,1,3,0),(2,0,2,1)
    one = [q for q in F if any(i==j for i,j in zip(q,vertex))]
    assert all(sum(i==j for i,j in zip(q,vertex))==1 for q in one)
    zero = [q for q in product(K,repeat=4) if all(i!=j for i,j in zip(q,vertex))]
    def shear(keys):
        return s.Matrix([[F.get(q[:p]+(vertex[p],)+q[p+1:],0) if q[p]==lift[p] else 0
                          for p in K] for q in keys])
    assert len(one)==28 and all(not F[q].has(b) for q in one)
    fixed = [q for q in one if not F[q].has(b) and s.Poly(s.fraction(F[q])[0],a,d).length()==1]
    candidates = []
    for sw in product((0,1),repeat=4):
        flags = tuple(v*i for i,v in zip((1,2,1,2),sw))
        if all(action(flags,q) in F for q in fixed):
            assert all(action(flags,q) in F for q in F)
            candidates.append(sw)
    assert candidates==[(0,0,0,0),(1,1,1,1)]
    if not equal:
        assert not any(shear([q for q in zero if q in F]))
        eql = exponent_matrix([q for q,v in F.items() if not v.free_symbols])
        full = exponent_matrix(F)
        pos = [next(q for q,v in F.items() if v==f) for f in (a,b,d)]
        rows = exponent_matrix(pos)
        for j in (1,2):
            diff = [k-l for k,l in zip(rows.tolist()[0],rows.tolist()[j])]
            assert fmpz_mat(eql.tolist()+[diff]).hnf().tolist()[:12]==eql.hnf().tolist()[:12]
        assert fmpz_mat(eql.tolist()+[rows.tolist()[0]]).hnf().tolist()[:13]==full.hnf().tolist()[:13]
        U = shear(zero)
        assert U*s.Matrix([-1,1,-1,1])==s.zeros(len(zero),1)
        const = s.Matrix([list(U.row(i)) for i in range(U.rows) if not any(v.free_symbols for v in U.row(i))])
        assert const.rank()==3
        return
    F = {q:v.subs(a,1) for q,v in F.items()}
    one_matrix = s.Matrix(exponent_matrix(one).tolist())
    diag,_,R = smith(one_matrix)
    rank = one_matrix.rank()
    assert [diag[j,j] for j in range(rank)]==[1]*11
    groups = defaultdict(list)
    for q in zero:
        char = s.Matrix(exponent_matrix([q]).tolist())*R
        groups[tuple(char[rank:])].append(q)
    used = []
    for char,ks in groups.items():
        U = shear(ks)
        f = s.Matrix([F.get(q,0) for q in ks])
        if not U.rank() and not any(f):
            continue
        image = [s.factor((v.T*f)[0]) for v in U.T.nullspace()]
        assert U.rank()==3
        assert s.Matrix([[v.expand().coeff(b,0),v.expand().coeff(b,1)] for v in image]).rank()==2
        assert s.gcd_list(list(char)) in (1,-1)
        used.append(char)
    assert len(used)==1


def cone():
    from verify_negative_component_point_graph import evaluation,z,x
    E = evaluation().subs({a:1,b:-1,c:0,d:1,z:0})
    T = s.Matrix([[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]])
    Q = s.diag(1,-1,-1,0)
    O = s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    mats = [E.row(i).jacobian(x).T*T for i in range(6)]
    assert all((M+M.T)/2==M[0,0]*Q for M in mats)
    assert any(M[0,0] for M in mats) and O.rank()==4
    # Find the same hyperplane of alternating forms as on the smooth quadrics.
    k = next(M for M in mats if M[0,0])
    assert all(sum((M-M[0,0]*k/k[0,0])[i,j]*O[i,j] for i in K for j in K)==0 for M in mats)


if __name__=='__main__':
    generators()
    flags(False)
    flags(True)
    cone()
    print('PASS: projective shifts, both filtered comparisons and the bare quadric cone limit')
