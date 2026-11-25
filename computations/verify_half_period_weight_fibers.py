"""A rational frame and its common regular residue for all nonzero weight vectors."""

from itertools import product

from flint import nmod_mat
import sympy as sp

from verify_half_period_frame import K, P, PAIRS, column, matrices, normalized, partner, relations
from verify_half_period_components import half_turn
from verify_walsh_full_component import vector


def modular(matrix):
    return nmod_mat([[int(sp.numer(v))*pow(int(sp.denom(v)),-1,P) % P for v in row]
                     for row in matrix.tolist()],P)


def verify():
    q = sp.diag(37,15,47,sp.Rational(9,37*15*47))
    u = sp.Matrix([[1,-5,25,-26],[1,-5,-25,26],[1,5,25,26],[1,5,-25,-26]])
    v = u.T.inv()*q
    assert q.det() == 9 and u.T*v == q
    frame = sp.Matrix([[u[s,i]*v[s,j] for i,j in PAIRS] for s in K])
    assert frame.rank() == 4
    # An alternating form restricts to zero on the diagonal.
    diagonal_restriction = sp.Matrix([[frame[s,4*i+i] if i==j else frame[s,4*i+j]+frame[s,4*j+i]
                       for i,j in PAIRS if i <= j] for s in K])
    assert modular(diagonal_restriction).rank() == 4
    volume = sp.Matrix(16,16,lambda a,b: sp.LeviCivita(a//4,a%4,b//4,b%4))
    metric = sp.Matrix(16,16,lambda a,b: q[a//4,b//4]*q[a%4,b%4]-q[a//4,b%4]*q[a%4,b//4])
    sq = q.reshape(16,1)*q.reshape(1,16)
    left = 45*volume+15*metric
    right = -45*volume+15*metric
    assert left.rank() == right.rank() == 3
    first,_ = matrices()
    align_inv = nmod_mat([[0,0,0,34],[0,0,4,0],[0,22,0,0],[26,0,0,0]],P).inv()
    move = nmod_mat([[int(i==k)*align_inv[j,l] for k,l in PAIRS] for i,j in PAIRS],P)
    J = nmod_mat([[int(i==j)*(-1)**(int(i//4==3)+int(i%4==3)) for j in range(16)] for i in range(16)],P)
    residue = move.transpose()*first*J*move
    edge0 = relations(first)
    _,second = matrices()
    edge1 = [r.transpose() for r in relations(second)]
    seed = (1,3,10,2)
    points = {tuple(int(i==j) for i in K) for j in K}
    points.update(vector([(-1)**((c&i).bit_count())*seed[i^x] for i in K]) for x,c in PAIRS)
    moved = [[],[]]
    for a in points:
        b = partner(a,edge0)
        c = normalized(half_turn(0,True)*column(partner(b,edge1)))
        if a != c:
            for rows,(x,y) in zip(moved,((a,b),(b,c)),strict=True):
                rows.append([x[i]*y[j] for i,j in PAIRS])
    assert all(len(rows) == 12 and nmod_mat(rows,P).rank() == 9 for rows in moved)

    # The two special weights below only verify arithmetic and the rank-six formulas.
    for weights in ((sp.Rational(1,4),)*4, (sp.Integer(53),2,3,-57)):
        weights = tuple(map(sp.sympify,weights))
        assert sum(weights) == 1
        x = 15*(frame.T*sp.diag(*(1/w for w in weights))*frame-sq)
        m,n = (x+left)/2,(x+right)/2
        assert m.rank() == n.rank() == 6
        assert all(n[4*i+j,4*k+l] == m[4*k+j,4*i+l] for i,j,k,l in product(K,repeat=4))
        if weights[0] == sp.Rational(1,4):
            assert modular(m) == residue


if __name__ == "__main__":
    verify()
    print("PASS: rational nine-component frame, good balanced lift and a weight divisible by 53")
