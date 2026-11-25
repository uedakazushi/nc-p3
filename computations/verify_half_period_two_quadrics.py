"""The order-four half-period component as a family of threefolds cut by two quadrics."""

from itertools import combinations

import sympy as sp

from verify_half_period_components import pencil_chart
from verify_walsh_full_component import P


def normalized_data(u, v):
    return (u + 1/u)/2, (u - 1/u)/2, (v + 1/v)/2, (v - 1/v)/2


def crossing_weights(h,k,j,l):
    swap = (0,2,1,3)
    first = sp.zeros(4)
    second = sp.zeros(4)
    for i,(a,b) in enumerate(zip((h,k,k,h),(j,l,l,j),strict=True)):
        first[i,swap[i] ^ 3] = a
        second[i,swap[i]] = b
    return first,second


def verify_transfer():
    u,v = sp.symbols("u v", nonzero=True)
    h,k,j,l = normalized_data(u,v)
    t = h*h-j*j
    a,b = crossing_weights(h,k,j,l)
    C = sp.Matrix(4,4,lambda i,j: int(i ^ j == 3))
    K = sp.diag(1,-1,-1,1)*C
    Z,I = sp.zeros(4),sp.eye(4)
    def blocks(a,b,c,d):
        return a.row_join(b).col_join(c.row_join(d))
    transfer = blocks(I,Z,-K,I)*blocks(Z,b.T.inv(),-b,Z)*blocks(I,Z,t*K,I)*blocks(Z,a.inv(),a.T,Z)
    transfer = transfer.applyfunc(sp.factor)
    assert (transfer**2 - sp.eye(8)).applyfunc(sp.factor) == sp.zeros(8)
    # On the two eigengraphs the second Schur transfer exchanges the symmetric offsets.
    D = sp.diag(j*h,-l*k,-l*k,j*h)/t
    A = sp.diag(h/j,k/l,k/l,h/j)
    R = -D*A*C
    for eps in (-1,1):
        graph = (R+eps*D).applyfunc(sp.factor)
        columns = I.col_join(graph)
        assert (transfer*columns - eps*columns).applyfunc(sp.factor) == sp.zeros(8,4)
        assert (a.T*graph.inv()*a-t*(R-eps*D)).applyfunc(sp.factor) == Z

    L = j*j+k*k
    eigen = (j*h,-j*h,l*k,-l*k)
    scales = ((j-h)**2,(j+h)**2,(l-k)**2,(l+k)**2)
    for i,m in combinations(range(4),2):
        weight = eigen[i]+eigen[m]
        assert sp.factor((weight+L)*scales[i]*scales[m] + t*t*(weight-L)) == 0


def verify_walsh_intersection():
    s,r,h,j = pencil_chart(1)
    hp = [next(h[i,k] for k in range(4) if h[i,k]) for i in range(4)]
    jp = [next(j[i,k] for k in range(4) if j[i,k]) for i in range(4)]
    h2,k2 = hp[0]*hp[3], hp[1]*hp[2]
    j2,l2 = jp[0]*jp[3], jp[1]*jp[2]
    ds,dt = s[1,2]-s[0,3],r[1,2]-r[0,3]
    assert ds*dt == h2-j2 == k2-l2 != 0
    assert h2-k2 == j2-l2 == 1
    assert (int(h2),int(k2),int(j2),int(l2)) == (44,43,37,36)
    assert (j2+k2)*(j2*h2-l2*k2)*(j2*l2-h2*k2) != 0
    from flint import nmod_mat
    I = nmod_mat([[int(i == j) for j in range(4)] for i in range(4)],P)
    K = nmod_mat([[(-1)**i.bit_count()*int(i ^ j == 3) for j in range(4)] for i in range(4)],P)
    on_graph = j.transpose().inv()*(h.transpose() + K*h.inv()*s*dt)
    assert (on_graph-I).rank() == (on_graph+I).rank() == 2
    a,b,c,d,e,f,g,v,w = sp.symbols("a b c d e f g v w")
    signed = sp.Matrix([[a*a,-b*b,-c*c],[-d*d,e*e,f*f],[-g*g,v*v,w*w]])
    rows,cols = signed*sp.ones(3,1),signed.T*sp.ones(3,1)
    assert sp.expand(c*c+v*v-a*a+e*e - (cols[1]-rows[0])) == 0
    assert sp.expand(f*f+g*g-a*a+e*e - (rows[1]-cols[0])) == 0


def verify_reconstruction():
    from flint import nmod, nmod_mat
    from verify_half_period_components import submatrix
    h,k,j,l = normalized_data(nmod(3,P), nmod(2,P))
    t = h*h-j*j
    d = (j*h/t,j*h/t,-l*k/t,-l*k/t)
    weights = (j*h,-j*h,l*k,-l*k)
    # In eigen-coordinates use rows (1,0),(0,1),(a,b),(c,e).
    first,second = crossing_weights(*(int(x) for x in (h,k,j,l)))
    first,second = [nmod_mat([[int(x) for x in row] for row in m.tolist()],P) for m in (first,second)]
    def diag(entries):
        return nmod_mat([[x*int(i == m) for m in range(4)] for i,x in enumerate(entries)],P)
    basis = nmod_mat([[1,1,0,0],[0,0,1,1],[0,0,1,-1],[1,-1,0,0]],P)
    C = nmod_mat([[int(i ^ m == 3) for m in range(4)] for i in range(4)],P)
    I = diag([1]*4)
    D = diag([d[0],d[2],d[2],d[0]])
    R = -D*diag([h/j,k/l,k/l,h/j])*C
    K = diag([1,-1,-1,1])*C
    for a in range(1,P):
        for e in range(1,P):
            U = nmod_mat([[1,0],[0,1],[a,1],[2,e]],P)
            minors = [U[i,0]*U[m,1]-U[m,0]*U[i,1] for i,m in combinations(range(4),2)]
            if sum(d[i]*d[m]*(weights[i]+weights[m]-j*j-k*k)*z*z
                   for z,(i,m) in zip(minors,combinations(range(4),2),strict=True)):
                continue
            U = basis*U
            gram = U.transpose()*D*U
            if not gram.det():
                continue
            S = R+D*(2*U*gram.inv()*U.transpose()*D-I)
            if not S.det() or not (S+K).det():
                continue
            T = first.transpose()*S.inv()*first
            assert S[1,2]-S[0,3] == 1 and T[1,2]-T[0,3] == t
            assert T+t*K == -second.transpose()*(S+K).inv()*second
            break
        else:
            continue
        break
    else:
        raise AssertionError("no plane with both Schur charts invertible")

    def indices(a,b):
        return [4*(2*a+i)+2*b+m for i in range(2) for m in range(2)]
    pos,neg = indices(0,0),indices(1,1)
    left,right = indices(0,1),indices(1,0)
    M = nmod_mat(16,16,P)
    for i,m in ((i,m) for i in range(4) for m in range(4)):
        M[pos[i],pos[m]],M[neg[i],neg[m]] = S[i,m],T[i,m]
        M[pos[i],neg[m]] = M[neg[m],pos[i]] = first[i,m]
    for labels in (left,right):
        M[labels[0],labels[3]],M[labels[3],labels[0]] = 1,-1
    for i in range(2):
        for m in range(2):
            entry = j if i == m else l
            M[left[3*i],right[3*m]] = entry
            M[right[3*m],left[3*i]] = -entry
    N = nmod_mat([[(-1)**(i//2)*M[4*i+b,4*k+n] for i in range(4) for n in range(4)]
                  for k in range(4) for b in range(4)],P)
    assert submatrix(N,pos,neg) == second
    for cut in (M,N):
        assert cut.rank() == 6
        assert (cut+cut.transpose()).rank() == 4
        assert (cut-cut.transpose()).rank() == 2


if __name__ == "__main__":
    verify_reconstruction()
    verify_transfer()
    verify_walsh_intersection()
    print("PASS: symplectic graph, Pluecker quadrics, smooth shared rank-fifteen point")
