"""Simultaneous rank completion, smoothing, length-four obstruction and positive AS certificates."""
from itertools import product

import sympy as s
from sympy.polys.matrices import DomainMatrix

from verify_klein_elliptic_completion import coeffs, four_relations, gram_twist
from verify_klein_period_four_center import K, artinian_quotient, commuting_space, paths

a, b, c, d, f, g, lam = s.symbols("a b c d f g l")
swap = s.Matrix([[0, 1], [1, 0]])
left = s.Matrix([[1, lam], [lam, 1]])
right = s.Matrix([[1, 1/lam], [1/lam, 1]])


def schur(m):
    return (left * m.inv() * right).applyfunc(s.factor)


def mixed_equations():
    A = s.Matrix([[a,b], [b,d]])
    B = s.Matrix([[b,f], [g,b]])
    T, U = schur(A), schur(B)
    Z = s.Matrix([[U[0,0], T[0,1]], [T[1,0], U[1,1]]])
    eq = s.factor(T[0,0] - U[0,1])
    assert T[0,0] == T[1,1] and U[0,1] == U[1,0]
    # The other three-by-three determinant follows everywhere on the localized hypersurface.
    second = b*Z.det() - Z[0,0]-Z[1,1]+Z[0,1]+Z[1,0]
    den = s.fraction(s.factor(second/eq))[1]
    assert not s.rem(s.Poly(lam*A.det()*B.det(),a,b,d,f,g,lam),
                     s.Poly(den,a,b,d,f,g,lam))
    return A, B, eq


def smoothing():
    delta, k = s.symbols("delta k")
    co, P, Q = coeffs(a,b,b+delta*k,d,1+delta,lam)
    def residue(expr):
        return s.factor(s.cancel(expr).subs(delta,0))
    slope = (1/lam-lam)/(2*(a*d-b*b))
    assert s.factor(residue(co[1]).diff(k)-slope) == 0
    assert s.factor(residue(co[1])+residue(co[3])+
                    (a+d-(lam+1/lam)*b)/(a*d-b*b)) == 0
    assert residue(delta*co[5]) != 0
    assert s.factor(residue(delta*co[5])+residue(delta*co[7])) == 0
    assert residue(delta**2*Q[0]) != 0
    assert s.factor(residue(delta*Q[1])) == 0


def opposite_flag_dimensions():
    A = s.Matrix([[a,b],[c,d]])
    B = A.T * swap
    def block(m,x,y):
        return m.row_join(x).col_join(y.row_join(y*m.inv()*x))
    blocks = block(A,right,left), block(B,lam*right,left/lam).T
    assert (schur(B)-schur(A).T*swap).applyfunc(s.factor) == s.zeros(2)
    def projection(h,transpose):
        if h >= 2:
            return s.eye(3).row_join(-s.ones(3,1))
        m = blocks[h].T if transpose else blocks[h]
        return (-(m[2:,:2]*m[:2,:2].inv())).row_join(s.eye(2)).applyfunc(s.factor)
    outer = [[projection(h,transpose) for h in K] for transpose in (True,False)]
    middle = [(h,[int(i//2==m) for i in K]) for h in (0,1) for m in (0,1)]
    middle += [(h,[1]*4) for h in (2,3)]
    dimensions = []
    for k in K:
        dims = [outer[0][h].rows*outer[1][h^k].rows for h in K]
        columns = []
        for i,j in product(K,repeat=2):
            for h,rel in middle:
                if i^j^h != k:
                    continue
                col = s.zeros(sum(dims),1)
                for u,r in enumerate(rel):
                    if r:
                        v = u^h
                        z = i^u
                        off = sum(dims[:z])
                        col[off:off+dims[z],:] += s.kronecker_product(
                            outer[0][z][:,i],outer[1][v^j][:,v])
                columns.append(col)
        rows = []
        for row in s.Matrix.hstack(*columns).tolist():
            fractions = [s.fraction(s.factor(v)) for v in row]
            den = s.lcm_list([denom for _,denom in fractions])
            rows.append([s.cancel(den*num/denom) for num,denom in fractions])
        rank = DomainMatrix.from_Matrix(s.Matrix(rows)).rank()
        dimensions.append(sum(dims)-rank)
    assert dimensions == [13,13,12,12]


def samples(A,B,eq):
    for coords in ((2,3,7,5,2,s.Rational(2,5)),(3,5,7,4,3,s.Rational(87,17))):
        sub = dict(zip((a,b,d,f,lam,g),map(s.sympify,coords)))
        assert eq.subs(sub)==0
        w,opp,n = gram_twist(A.subs(sub),B.subs(sub),s.S.One,sub[lam])
        assert [m.rank() for m in (w,opp,n)] == [6,16,6]
        rel = four_relations(w,n)
        basis,image = paths(rel,5)
        assert all([sum(len(basis[p,t,h]) for h in K) for t in range(6)]
                   == [1,4,10,20,35,56] for p in K)
        assert [commuting_space(basis,image,h,0).cols for h in K] == [3,2,0,0]
        artinian_quotient(rel,basis,image)


if __name__ == "__main__":
    A,B,eq = mixed_equations()
    smoothing()
    opposite_flag_dimensions()
    samples(A,B,eq)
    print("PASS: mixed resolution, smoothing residues, all opposite flags excluded and two central cuts")
