"""Verify the forced splitting, the compressed three-letter block, and the exponential bound."""
from itertools import combinations_with_replacement,product
from math import comb

import sympy as sp
from sympy.polys.matrices import DomainMatrix


def main():
    possible = [a for a in combinations_with_replacement(range(1,8),6)
                if sum(a)==12 and sum(max(d-3,0) for d in a)>=4]
    assert possible == [(1,1,1,1,1,7)]
    splitting = possible[0]
    assert [sum(max(d-i-1,0) for d in splitting) for i in range(7)] == [6,5,4,3,2,1,0]
    # A generic four-dimensional octavic system. Its seven shifts give both
    # coefficients of a binary diagonal on contraction with degree-five tests.
    S = sp.Matrix(9,4,lambda i,j: ((i+2)**(j+2)+13*i*j+5*j*j)%37)
    basis = list(combinations_with_replacement(range(3),2))
    cols = []
    for (i,j),k in product(basis,range(4)):
        M = sp.zeros(3,7)
        for row,h in {(i,j),(j,i)}:
            for b in range(7):
                M[row,b] = S[h+b,k]
        cols.append(M.reshape(21,1))
    M = sp.Matrix.hstack(*cols)
    assert 24-DomainMatrix.from_Matrix(M).rank() >= 3
    for vector in M.nullspace():
        triple = sp.MutableDenseNDimArray.zeros(3,4,3)
        for a,((i,j),k) in zip(vector,product(basis,range(4))):
            triple[i,k,j] = triple[j,k,i] = a
        for side,delta,index,test in product(range(2),range(2),range(3),range(6)):
            assert sum((triple[index,k,h] if side==0 else triple[h,k,index])*S[h+delta+test,k]
                       for k in range(4) for h in range(3)) == 0
    assert 4*3**7 > comb(32,3)
    print("Unique splitting (1,1,1,1,1,7), both diagonal contractions and the length-29 Hilbert obstruction verified.")


if __name__ == '__main__':
    main()
