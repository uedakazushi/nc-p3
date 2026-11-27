"""Thin quotient under four-letter relabeling and the final invertible-corner fiber."""
from itertools import product

import sympy as s
from flint import fmpz_mat
from sympy.matrices.normalforms import hermite_normal_form

from verify_double_exchange_fourier import twist
from verify_exchange_singular_corners import index_row
from verify_exchange_mixed_configurations import incidence_signatures, maps, modular_points, partner
from verify_klein_elliptic_completion import four_relations, gram_twist
from verify_klein_period_four_center import artinian_quotient, commuting_space, paths
from verify_negative_component_quotient import K, lattice, relabel, thin, words
from verify_parallel_invertible_rectangles import smith

z,lam = s.symbols('z l')


def tensor(A,B,z,lam):
    w,_,_ = gram_twist(A,B,z,lam)
    return {q:s.factor(w[4*q[0]+q[1],4*q[2]+q[3]]) for q in words}


def thin_group():
    A = s.Matrix([[2,3],[5,7]])
    F = tensor(A,A+3*s.eye(2),z,lam)
    ex = {q:[int(F[q].as_powers_dict().get(v,0)) for v in (z,lam)] for q in thin}
    ker = lattice(thin)
    weights = s.Matrix((ker*fmpz_mat([ex[q] for q in thin])).tolist())
    assert hermite_normal_form(weights.T)==s.eye(2)
    actions = [set(),set()]
    for shift,shear,i in product(K,range(2),K):
        for j,k in product((i&2,(i&2)^1),repeat=2):
            g = shift,shear,(i,j,k,i^j^k)
            v = s.Matrix((ker*fmpz_mat([ex[relabel(q,*g)] for q in thin])).tolist())
            change = weights.gauss_jordan_solve(v)[0]
            assert weights*change==v
            actions[shift%2].add(tuple(change))
    assert actions[0]=={(1,0,0,1),(1,0,-1,-1),(-1,0,0,-1),(-1,0,1,1)}
    assert actions[1]=={(1,2,-1,-1),(1,2,0,-1),(-1,-2,1,1),(-1,-2,0,1)}


def special_actions():
    a,b,c,d,e,f,g,h = u = s.symbols('a b c d e f g h')
    A,B = s.Matrix(2,2,u[:4]),s.Matrix(2,2,u[4:])
    F = tensor(A,B,-s.S.One,s.S.One)
    E = s.Matrix([index_row(q) for q in thin])
    H,L,R = smith(E)
    rank = E.rank()
    for shift,tw,off,up in [(1,0,(0,0,1,1),(g,a,-d,-f,c,e,-h,-b)),
                            (0,1,(0,0,1,1),(b,a,-d,-c,g,-h,e,-f)),
                            (2,0,(0,0,0,0),(a,c,b,d,e,g,f,h)),
                            (0,0,(1,0,1,0),(h,-f,-g,e,d,-b,-c,a))]:
        G = {q:F[relabel(q,shift,0,off)]*twist(q)**tw for q in words}
        rhs = L*s.Matrix([2*int(G[q]!=F[q]) for q in thin])
        assert all(v%4==0 for v in rhs[rank:])
        phase = R*s.Matrix([rhs[j]/H[j,j] if j<rank else 0 for j in range(16)])
        image = {q:v*s.I**(int(sum(phase[4*p+i] for p,i in enumerate(q)))%4) for q,v in G.items()}
        assert all(image[q]==F[q] for q in thin)
        top = [image[i,i,j,j] for i,j in product(range(2),repeat=2)]
        top += [image[j,j^1,i,i^1] for i,j in product(range(2),repeat=2)]
        assert tuple(top) in (up,tuple(-v for v in up))
    X = s.Matrix([[1,1],[1,-1]])
    Y = s.Matrix([[1,1],[-1,1]])
    T = Y*A.adjugate()*X
    assert T==s.Matrix([[a+d-b-c,-a+d+b-c],[a-d+b-c,-a-d-b-c]])


def positive_fibers():
    examples = [([[2,3],[5,7]],[[-7,3],[5,-2]]),([[2,3],[3,7]],[[5,3],[3,1]])]
    for A,B in examples:
        w,o,n = gram_twist(s.Matrix(A),s.Matrix(B),-s.S.One,s.S.One)
        assert [m.rank() for m in (w,o,n)]==[6,10,6]
        rel = four_relations(w,n)
        basis,image = paths(rel,5)
        assert all([sum(len(basis[p,t,h]) for h in K) for t in range(6)]
                   ==[1,4,10,20,35,56] for p in K)
        assert [commuting_space(basis,image,h,0).cols for h in K]==[3,2,0,0]
        artinian_quotient(rel,basis,image)


def special_maps():
    p = 109
    w,_,n = gram_twist(s.Matrix([[2,3],[5,7]]),s.Matrix([[-7,3],[5,-2]]),-s.S.One,s.S.One)
    pts = [modular_points(w,p)]
    for mat,tr in ((w,False),(n,True),(w.T,False)):
        pts.append(partner(pts[-1],mat,tr,p))
    assert partner(pts[-1],n.T,True,p)==pts[0]
    flip = [4*j+i for i,j in product(K,repeat=2)]
    for points,mat in zip(pts[1:],(n[flip,:],w.T,n.T[flip,:])):
        assert set(modular_points(mat,p))==set(points)
    tagged = []
    for points in pts:
        hist,sigs = incidence_signatures(points,p)
        tagged.append((points,[(tuple(sorted(hist.items())),sig) for sig in sigs]))
    assert [[maps(x,y,p) for y in tagged] for x in tagged]==[
        [8,0,8,0],[0,8,0,0],[8,0,8,0],[0,0,0,8]]


if __name__=='__main__':
    thin_group()
    special_actions()
    positive_fibers()
    special_maps()
    print('PASS: both marked thin action groups, two last-fiber reductions and central certificates')
