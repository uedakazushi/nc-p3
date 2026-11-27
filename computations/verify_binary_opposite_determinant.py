"""The leading opposite determinant along both symmetry divisors and its derivative."""
from itertools import permutations

import numpy as np
from flint import nmod_mat

from verify_binary_component_crossing import (
    PRIME, character_matrices, mod_kernel, reduced, sample, tangent_data,
)


def matrix(array):
    return nmod_mat(array.tolist(),PRIME)


def multiply(*factors):
    result = factors[0]
    for factor in factors[1:]:
        result = result@factor%PRIME
    return result


def pivots(array):
    echelon, rank = matrix(array).rref()
    return [next(j for j in range(echelon.ncols()) if echelon[i,j]) for i in range(rank)]


def parity(sign):
    columns = []
    for i in range(4):
        for j in range(i if sign==1 else i+1,4):
            v = np.zeros(16,dtype=int)
            v[4*i+j],v[4*j+i] = 1,sign
            columns.append(v)
    return np.array(columns).T


def opposite(array):
    return array.reshape(4,4,4,4).transpose(0,2,1,3).reshape(16,16)


def main():
    w, integral, n, inverse = sample()
    symmetric, mixed, forms = tangent_data(integral,n,inverse)
    A, B = mixed[1,-1], mixed[-1,1]
    def cross(X,Y):
        return (multiply(X,inverse,Y)+multiply(Y,inverse,X))%PRIME

    def solve(character, source):
        basis = character_matrices(*character)
        idx = np.triu_indices(10,0 if character[0]*character[1]==1 else 1)
        jac = multiply(n.T,basis,n)[:,idx[0],idx[1]].T
        rhs = multiply(n.T,source,n)[idx]
        cols = pivots(jac)
        rows = pivots(jac[:,cols].T)
        inv = np.array(matrix(jac[np.ix_(rows,cols)]).inv().tolist(),dtype=int)
        values = multiply(inv,rhs[rows])
        assert np.all((jac[:,cols]@values-rhs)%PRIME==0)
        return (values@basis[cols].reshape(-1,256)).reshape(16,16)%PRIME

    Z = solve((-1,-1),cross(A,B))
    sy, sk = parity(1), parity(-1)
    def scaled_cut(base,first,second,both):
        return np.block([[multiply(sy.T,opposite(base),sy),multiply(sy.T,opposite(first),sk)],
                         [multiply(sk.T,opposite(second),sy),multiply(sk.T,opposite(both),sk)]])
    leading = scaled_cut(reduced(w),A,B,Z)
    O = leading[:10,:10]
    assert matrix(mod_kernel(O.T).T@leading[:10,10:]%PRIME).rank() == 1
    assert matrix(leading[10:,:10]@mod_kernel(O)%PRIME).rank() == 2
    left, right = mod_kernel(leading.T), mod_kernel(leading)
    assert left.shape[1] == right.shape[1] == 1

    # A symmetric direction tangent to the common smooth intersection.
    X = (mod_kernel(forms)[:,0]@symmetric.reshape(45,256)).reshape(16,16)%PRIME
    first = solve((1,-1),cross(X,A))
    second = solve((-1,1),cross(X,B))
    source = cross(X,Z)+cross(A,second)+cross(B,first)
    for U,V,W in permutations((X,A,B)):
        source -= multiply(U,inverse,V,inverse,W)
    third = solve((-1,-1),source)
    derivative = scaled_cut(X,first,second,third)
    assert multiply(left.T,derivative,right)[0,0] != 0
    print("Ranks 1 and 2, scaled opposite rank 15, and a transverse derivative of its determinant verified.")


if __name__ == '__main__':
    main()
