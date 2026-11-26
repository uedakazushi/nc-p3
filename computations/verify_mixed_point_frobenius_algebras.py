"""Check the mixed construction over Q(i) on two nonreduced Frobenius algebras."""
from itertools import product

from sympy import I, Matrix, eye, zeros


def check(local_type):
    basis = [eye(4)[:, j] for j in range(4)]
    unit, socle = basis[0], basis[3]

    def multiply_k(a, b):
        result = zeros(4, 1)
        for i, j in product(range(4), repeat=2):
            if local_type == "curvilinear":
                k = i + j
            elif not i or not j:
                k = i + j
            else:
                k = 3 if {i, j} == {1, 2} else 4
            if k < 4:
                result[k] += a[i]*b[j]
        return result

    def tau(a):
        return a[0] + a[3]
    U = Matrix.hstack(*basis[1:]) - unit*Matrix([[0, 0, 1]])
    gram = Matrix(3, 3, lambda i,j: -tau(multiply_k(U[:,i], U[:,j])))
    assert gram.det() != 0

    def multiply(a, b):
        r, s = a[0], b[0]
        x, y = a[1:4, :], b[1:4, :]
        return Matrix([r*s-x.dot(y)]).col_join(r*y+s*x-x.cross(y)).col_join(
            multiply_k(a[4:,:], b[4:,:]))

    def span(columns):
        return Matrix.hstack(*columns).columnspace()
    def products(left, right):
        return span([multiply(a,b) for a,b in product(left,right)])

    even = [unit.col_join(zeros(4,1))]
    even.extend(basis[i+1].col_join(U[:,i]) for i in range(3))
    odd = [zeros(4,1).col_join(unit)]
    odd.extend(basis[i+1].col_join(U*gram.inv()[:,i]) for i in range(3))
    hyper = products(even,odd)
    assert len(hyper) == len(products(odd,even)) == 7
    assert all(a[0] == tau(a[4:,:]) for a in hyper + products(odd,even))
    assert len(products(hyper,even)) == len(products(hyper,odd)) == 8
    quotient = lambda a: a[1:4,:].col_join((a[4:,:]-tau(a[4:,:])*unit)[1:,:])
    assert len(span(list(map(quotient,hyper)))) == 6

    # A self-annihilating length-two ideal. Each pencil has split roots over Q(i).
    i = 2 if local_type == "curvilinear" else 1
    ideal = [basis[i],socle]
    assert not any(any(multiply_k(a,b)) for a,b in product(ideal,ideal))
    def odd_lift(k):
        return Matrix([0]).col_join(gram*(k-tau(k)*unit)[1:,:]).col_join(k)
    even_pair = [even[0],even[i]]
    odd_pair = list(map(odd_lift,ideal))
    assert odd_lift(socle)[1:4,:].dot(odd_lift(socle)[1:4,:]) != 0
    full = [eye(8)[:,j] for j in range(8)]
    for side, partner, nxt in ((even_pair,odd_pair,odd),(odd_pair,even_pair,even)):
        from sympy import solve, symbols
        t = symbols("t")
        pencil = side[0]+t*side[1]
        roots = [pencil.subs(t,r) for r in solve(pencil[:4,:].dot(pencil[:4,:]),t)]
        assert roots
        for root in roots:
            coefficients = Matrix.hstack(*(multiply(root,r) for r in partner)).nullspace()
            assert len(coefficients) == 1
            paired = Matrix.hstack(*partner)*coefficients[0]
            assert not any(multiply(root,paired))
            assert [len(products([root],space)) for space in (nxt,hyper,full)] == [3,4,4]
        assert [len(products(side,space)) for space in (nxt,hyper,full)] == [5,6,6]

    net = []
    for j in range(3):
        n = gram.inv()[:,j]
        B = Matrix.hstack(*[(multiply_k(U*n,U[:,i])-tau(multiply_k(U*n,U[:,i]))*unit)[1:,:]
                            for i in range(3)])
        assert gram*B == B.T*gram
        net.append(B)
    for i,j in product(range(3),repeat=2):
        a,b = gram.inv()[:,i],gram.inv()[:,j]
        assert net[i]*net[j]-net[j]*net[i] == (a*b.T-b*a.T)*gram


if __name__ == "__main__":
    check("curvilinear")
    check("embedding-dimension-two")
    print("PASS: two nonreduced quartic algebras, cyclic products, self-annihilating ideals and mixed commutators")
