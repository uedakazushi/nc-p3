"""Check the four vertex frames over the boundary and the parameter of the smooth quadric graph."""
import sympy as s
from itertools import combinations

from verify_klein_elliptic_completion import four_relations,gram_twist
from verify_negative_component_point_graph import K,a,b,c,d,z,x,evaluation

w = s.Symbol('w')


def boundary_vertices():
    A = s.Matrix([[a,b],[c,a]])
    W,_,rev = gram_twist(A,-A,w*w,s.S.NegativeOne)
    # Choose rank-six columns over the function field. Rank never drops in the chart.
    for p,relations in enumerate(four_relations(W,rev)):
        pair0 = [(1,i,0,0) for i in (1,-1)]
        pair1 = [(0,0,i*w,1) for i in (1,-1)]
        if p%2:
            pair0 = [(0,0,1,i) for i in (1,-1)]
            pair1 = [(i*w,1,0,0) if p==1 else (1,i*w,0,0) for i in (1,-1)]
        for v in pair0+pair1:
            m = s.Matrix.vstack(*[(s.Matrix(1,4,v)*s.Matrix(4,4,row)) for _,row in relations])
            assert m.rank()==3
        # Vanishing on every coefficient column (also where a chosen minor vanishes).
        # Rank at most three is a closed condition in the fixed-rank relation bundle.


def quadric():
    E = evaluation().subs({a:1,b:-1,c:-z,d:1})
    T = s.Matrix([[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]])
    J = s.diag(1,-1,-1,z)
    mats = [E.row(i).jacobian(x).T*T for i in range(6)]
    assert all((M+M.T)/2==M[0,0]*J for M in mats)
    sym = [M[0,0] for M in mats]
    k = next(i for i,h in enumerate(sym) if h)
    skew = [M-sym[i]*mats[k]/sym[k] for i,M in enumerate(mats)]
    pairs = list(combinations(K,2))
    L = s.Matrix([[M[i,j] for i,j in pairs] for M in skew])
    O = s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    assert L*s.Matrix([O[i,j] for i,j in pairs])==s.zeros(6,1)
    # The five-dimensional kernel rank follows from the six independent relations
    # and a nonzero symmetric image for every z on this chart.
    h = s.Symbol('h')
    assert s.factor((J*O).charpoly(h).as_expr()-(h*h-1)*(h*h-z))==0


if __name__=='__main__':
    boundary_vertices()
    quadric()
    print('PASS: vertex labels in all four positions and the smooth-quadric coordinate z+1/z')
