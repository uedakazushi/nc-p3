"""The eight fixed points, the two rulings, and the four rank-one bilinear forms."""

from itertools import combinations_with_replacement, product

from flint import nmod_mat
from sympy import LeviCivita

from verify_half_period_components import half_turn
from verify_walsh_full_component import K, P, matrices, vector

PAIRS = list(product(K,repeat=2))
SYM = list(combinations_with_replacement(K,2))


def column(v):
    return nmod_mat([[int(a)] for a in v],P)


def normalized(v):
    return vector([int(a) for a in v.entries()])


def relations(m):
    rows,rank = m.transpose().rref()
    assert rank == 6
    return [nmod_mat(4,4,rows.entries()[16*i:16*(i+1)],P) for i in range(rank)]


def partner(a, side):
    row = column(a).transpose()
    kernel,dim = nmod_mat([list((row*r).entries()) for r in side],P).nullspace()
    assert dim == 1
    return vector([int(kernel[i,0]) for i in K])


def verify():
    first,second = matrices()
    edge0,edge1 = relations(first),[r.transpose() for r in relations(second)]
    units = [tuple(int(i == j) for i in K) for j in K]
    seed = (1,3,10,2)
    points = set(units) | {vector([(-1)**((c&i).bit_count())*seed[i^x] for i in K]) for x,c in PAIRS}
    odd = {p:partner(p,edge0) for p in points}
    even = {p:partner(odd[p],edge1) for p in points}
    for x in K:
        perm = {p:normalized(half_turn(x,True)*column(even[p])) for p in points}
        assert set(perm.values()) == points
        fixed = {p for p in points if p == perm[p]}
        fixed_square = {p for p in points if p == perm[perm[p]]}
        assert (len(fixed),len(fixed_square)) == ((0,4) if x else (8,20))
        assert all(perm[perm[perm[perm[p]]]] == p for p in points)
        moved_rows = [[p[i]*odd[p][j] for i,j in PAIRS] for p in points-fixed]
        assert nmod_mat(moved_rows,P).rank() == (10 if x else 9)
    turn = {p:normalized(half_turn(0,True)*column(even[p])) for p in points}
    fixed = {p for p in points if p == turn[p]}
    frame = nmod_mat([list(odd[p]) for p in units],P).transpose()
    p = next(p for p in sorted(fixed) if all(p))
    b = frame.inv()*column(odd[p])
    align = nmod_mat([[int(i==j)*p[i]/b[i,0] for j in K] for i in K],P)*frame.inv()
    odd = {p:normalized(align*column(b)) for p,b in odd.items()}
    # The polarization of the common smooth quadric vanishes on both sets of secants.
    def quad(s,t):
        return [s[i]*t[j]+int(i!=j)*s[j]*t[i] for i,j in SYM]
    lines = []
    triples = [[],[]]
    for p in sorted(points):
        b,c,d = odd[p], turn[p], odd[turn[p]]
        if p in fixed:
            assert p == b == c
        else:
            for u,v in ((p,b),(b,c)):
                lines.extend((quad(u,u),quad(u,v),quad(v,v)))
        for rows,(u,v,w) in zip(triples,((p,b,c),(b,c,d)),strict=True):
            rows.append([u[i]*v[j]*w[k] for i,j,k in product(K,repeat=3)])
    assert all(nmod_mat(rows,P).rank() == 20 for rows in triples)
    net,dim = nmod_mat(lines,P).nullspace()
    assert dim == 1
    q = nmod_mat([[net[SYM.index(tuple(sorted((i,j)))),0] for j in K] for i in K],P)
    assert q.det() == 9
    assert all((column(p).transpose()*q*column(p))[0,0] for p in fixed)
    J = nmod_mat([[int(i==j)*(-1)**(int(i//4==3)+int(i%4==3)) for j in range(16)] for i in range(16)],P)
    gi = align.inv()
    move = nmod_mat([[int(i==k)*gi[j,l] for k,l in PAIRS] for i,j in PAIRS],P)
    m = move.transpose()*first*J*move
    n = nmod_mat([[m[4*k+j,4*i+l] for k,l in PAIRS] for i,j in PAIRS],P)
    volume = nmod_mat([[int(LeviCivita(i,j,k,l)) for k,l in PAIRS] for i,j in PAIRS],P)
    metric = nmod_mat([[q[i,k]*q[j,l]-q[i,l]*q[j,k] for k,l in PAIRS] for i,j in PAIRS],P)
    square = nmod_mat([[q[i,j]*q[k,l] for k,l in PAIRS] for i,j in PAIRS],P)
    h = (m-n)[1,11]
    assert m-n == h*volume and h == 45
    scalar = h/3
    X = m-(m-n+scalar*metric)*pow(2,-1,P)
    Y = X+scalar*square*pow(2,-1,P)
    assert X.rank() == 3 and Y.rank() == 4
    assert all(Y[4*i+j,4*k+l] == Y[4*k+j,4*i+l] for i,j,k,l in product(K,repeat=4))
    rows,nrows = Y.rref()
    forms = [nmod_mat(4,4,rows.entries()[16*i:16*(i+1)],P) for i in range(nrows)]
    def evaluate(a):
        return nmod_mat([list((column(a).transpose()*f).entries()) for f in forms],P)
    basis = evaluate((1,1,1,1))
    assert basis.det() == 3
    operator = basis.inv()*evaluate((7,4,2,9))
    assert sorted((int(p[0]),e) for p,e in operator.charpoly().factor()[1]) == [(3,1),(4,1),(15,1),(52,1)]


if __name__ == "__main__":
    verify()
    print("PASS: half-period cycles, eight fixed points, smooth ruling quadric, four simple eigenspaces")
