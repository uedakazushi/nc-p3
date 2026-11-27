"""Exact transport of an obstruction through both orientations of a common linear factor."""
from math import comb, factorial

import sympy as sp

from verify_binary_secant_obstructions import rank


def main():
    T = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    G = sp.Matrix(6, 6, lambda i,j: T[i+j])
    Y = [2, 3]
    for j in range(7):
        Y.append(Y[-1]+T[j+4]+2*T[j+1])
    L = G.inv()*sp.Matrix(6,4,lambda i,j:Y[i+j])
    x = sp.symbols("x")
    d, linear, cubic = x, x-1, x**3+2

    def coeff(form, length):
        return sp.Matrix([sp.expand(form).coeff(x,i) for i in range(length)])

    def hankel(matrix):
        return all(matrix[i,j+1] == matrix[i+1,j]
                   for i in range(matrix.rows-1) for j in range(matrix.cols-1))

    columns = [coeff(linear*x**j,6) for j in range(5)] + [coeff(cubic*x**j,6) for j in range(3)]
    images = [coeff(d*x**j,6) for j in range(5)] + [L*coeff(d*x**j,4) for j in range(3)]
    A, B = sp.Matrix.hstack(*columns), sp.Matrix.hstack(*images)
    pivots = A.rref()[1]
    transport = B[:,pivots]*A[:,pivots].inv()
    assert transport*A == B and transport.det() == 7
    G2 = transport.T*G*transport
    J = sp.Matrix.hstack(*(coeff(linear*x**j,4) for j in range(3)),coeff(cubic,4))
    extra = sp.eye(6)[:,:3].row_join(L[:,0])
    H2 = G2*transport.inv()*extra*J.inv()
    assert hankel(G2) and hankel(H2)

    def transvectant(left, right, order):
        m, n = len(left)-1, len(right)-1
        def falling(i,j):
            return factorial(i)//factorial(i-j) if 0 <= j <= i else 0
        result = [0]*(m+n-2*order+1)
        for i,c in enumerate(left):
            for j,h in enumerate(right):
                for k in range(order+1):
                    value = ((-1)**k*comb(order,k)*falling(m-i,order-k)*falling(i,k)
                             *falling(n-j,k)*falling(j,order-k))
                    if value:
                        result[i+j-order] += value*c*h
        return result

    def absolute_invariant(moments):
        form = [comb(10,i)*v for i,v in enumerate(moments)]
        quartic = transvectant(form,form,8)
        return sp.Rational(transvectant(quartic,quartic,4)[0],
                           transvectant(form,form,10)[0]**2)

    second_moments = [G2[min(i,5),max(i-5,0)] for i in range(11)]
    assert absolute_invariant(T) != absolute_invariant(second_moments)
    a,b = sp.symbols("a b")
    other = (3, 1, 4, 1, 5, 9, 2, 6, 5)
    pencil = sp.Matrix(5,7,lambda i,j:T[i+j]).nullspace()
    octavic = [a*u+b*v for u,v in zip(Y,other)]
    rows = sp.Matrix(2,3,lambda i,j:sum(pencil[i][k]*octavic[j+k] for k in range(7)))
    ell = rows.row(0).cross(rows.row(1)).T
    polar = sp.Matrix(7,5,lambda i,j:T[i+j])
    rhs = sp.Matrix(7,3,lambda i,j:octavic[i+j])*ell
    quartic = polar[:5,:].inv()*rhs[:5,:]
    assert (polar*quartic-rhs).applyfunc(sp.expand).is_zero_matrix
    def form(v):
        return sum(c*x**i for i,c in enumerate(v))
    resultant = sp.Poly(sp.resultant(form(ell), form(quartic), x), a,b)
    assert resultant.is_homogeneous and resultant.total_degree() == 14
    assert resultant.eval(b,0).is_zero and not resultant.diff(b).eval(b,0).is_zero
    assert rank(rows.subs({a:1,b:0})) == 2
    print("Both binary Gram presentations, transport determinant 7, and a transverse degree-fourteen resultant verified.")


if __name__ == '__main__':
    main()
