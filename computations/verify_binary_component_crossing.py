"""Exact tangent ranks and good-prime mixed Schur forms at a stable nonregular point."""
from itertools import product
from math import lcm

import numpy as np
import sympy as sp
from flint import fmpz_mat, nmod_mat


PRIME = 1009


def reduced(matrix):
    return np.array([int(v.p)*pow(int(v.q),-1,PRIME)%PRIME for v in matrix]).reshape(matrix.shape)


def mod_kernel(matrix):
    basis, dimension = nmod_mat(matrix.tolist(),PRIME).nullspace()
    return np.array(basis.tolist(),dtype=int)[:,:dimension]


def character_matrices(even, odd):
    matrices = []
    for i,j,k,l in product(range(4),repeat=4):
        if i>k or j>l or (i==k and even==-1) or (j==l and odd==-1):
            continue
        m = np.zeros((16,16),dtype=int)
        for u,v,y,z,c in ((i,j,k,l,1),(k,j,i,l,even),(i,l,k,j,odd),(k,l,i,j,even*odd)):
            m[4*u+v,4*y+z] = c
        matrices.append(m)
    return np.array(matrices)


def main():
    T = (2,3,5,7,11,13,17,19,23,29,31)
    Y = (3,1,4,1,5,9,2,6,5)
    G = sp.Matrix(6,6,lambda i,j:T[i+j])
    L = G.inv()*sp.Matrix(6,4,lambda i,j:Y[i+j])
    shifts = [sp.Matrix(6,4,lambda i,j:int(i==j+k)) for k in range(3)]
    relation = sp.Matrix.hstack(*shifts,L).T
    w = relation*G*relation.T

    sextics = sp.Matrix(5,7,lambda i,j:T[i+j]).nullspace()
    rows = sp.Matrix(2,3,lambda i,j:sum(sextics[i][k]*Y[j+k] for k in range(7)))
    ell = rows.row(0).cross(rows.row(1)).T
    Yell = sp.Matrix(7,3,lambda i,j:Y[i+j])*ell
    c = sp.Matrix(7,5,lambda i,j:T[i+j])[:5,:].inv()*Yell[:5,:]
    x = sp.symbols("x")
    assert ell[-1] and c[-1]
    assert sp.resultant(sum(u*x**i for i,u in enumerate(ell)),sum(u*x**i for i,u in enumerate(c)),x)

    N = sp.Matrix.hstack(*(v*lcm(*(int(t.q) for t in v)) for v in relation.T.nullspace()))
    integral = np.array([int(t) for t in N],dtype=object).reshape(N.shape)
    n = reduced(N)
    pivots = list(relation.T.rref()[1])
    inverse = np.zeros((16,16),dtype=int)
    inverse[np.ix_(pivots,pivots)] = reduced(w[pivots,pivots].inv())

    lifts, obstructions = {}, {}
    expected = {(1,1):45,(1,-1):16,(-1,1):16,(-1,-1):0}
    for (even,odd),dimension in expected.items():
        basis = character_matrices(even,odd)
        index = np.triu_indices(10,0 if even*odd==1 else 1)
        jacobian = (integral.T@basis@integral)[:,index[0],index[1]].T
        assert basis.shape[0]-fmpz_mat(jacobian.tolist()).rank() == dimension
        jacobian = np.asarray(jacobian%PRIME,dtype=int)
        tangent = mod_kernel(jacobian)
        assert tangent.shape[1] == dimension
        lifts[even,odd] = (tangent.T@basis.reshape(-1,256)).reshape(-1,16,16)%PRIME
        obstructions[even,odd] = mod_kernel(jacobian.T),index

    symmetric = lifts[1,1]
    forms = []
    for character in ((1,-1),(-1,1)):
        functional,index = obstructions[character]
        assert functional.shape[1] == 1
        pairings = []
        for mixed in lifts[character]:
            cross = ((symmetric@inverse%PRIME)@mixed+(mixed@inverse%PRIME)@symmetric)%PRIME
            normal = ((n.T@cross%PRIME)@n)%PRIME
            pairings.append(normal[:,index[0],index[1]]@functional[:,0]%PRIME)
        assert nmod_mat(np.array(pairings).tolist(),PRIME).rank() == 1
        forms.extend(pairings)
    assert nmod_mat(np.array(forms).tolist(),PRIME).rank() == 1
    print("Finite sample, four rational Jacobians, and nonzero mixed Hessians in both characters verified.")


if __name__ == '__main__':
    main()
