"""Verify the displacement identities for the binary obstruction, including the rank-one locus."""
import sympy as sp

from verify_binary_secant_obstructions import rank


def main():
    T = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    G = sp.Matrix(6, 6, lambda i,j: T[i+j])
    sextics = sp.Matrix(5, 7, lambda i,j: T[i+j])
    f = next(v / v[6] for v in sextics.nullspace() if v[6] != 0)
    e = G.inv()[:, 5].col_join(sp.zeros(1, 1))

    def coordinates(Y):
        H = sp.Matrix(6, 4, lambda i,j: Y[i+j])
        L = G.inv() * H
        rows = sp.Matrix(2, 3, lambda i,j: sum((f,e)[i][k] * Y[k+j] for k in range(7)))
        for j in range(3):
            defect = L[:,j+1].col_join(sp.zeros(1,1)) - sp.zeros(1,1).col_join(L[:,j])
            assert defect == rows[0,j]*e - rows[1,j]*f
        K = L.T * G * L
        for i in range(3):
            for j in range(3):
                assert K[i+1,j]-K[i,j+1] == rows[0,i]*rows[1,j]-rows[1,i]*rows[0,j]
        return H,L,rows

    basis = list(sp.eye(9).columnspace())
    samples = [coordinates(Y) for Y in basis]
    linear = sp.Matrix.hstack(*(a.reshape(6,1) for _,_,a in samples))
    assert rank(linear) == 6
    for expected in (1,2):
        target = sp.Matrix([1,0,0,0,expected-1,0])
        solution, freedom = linear.gauss_jordan_solve(target)
        Y = solution.subs({x: 0 for x in freedom})
        _, L, rows = coordinates(Y)
        assert rank(rows) == expected
        system = []
        for H, _, _ in samples:
            mixed = L.T * H
            system.append(sp.Matrix([mixed[i,j]-mixed[j,i] for i in range(4) for j in range(i+1,4)]))
        assert 9-rank(sp.Matrix.hstack(*system)) == (6 if expected==1 else 4)
    print("Both apolar coordinates, all displacement minors, and cubic kernel dimensions verified.")


if __name__ == '__main__':
    main()
