#!/usr/bin/env python3
"""Exact finite certificates used in the rank-two pencil exclusion proof."""
from fractions import Fraction as F
from itertools import product
from sympy import Matrix, eye, zeros, diag, kronecker_product
from frobenius_tensor import (flatten, reconstructed_windows,
    orthogonal_windows, nullspace_sparse)
from verify_quartic_obstruction import insertion_columns, exact_rank


def jordan(sizes, eigenvalue):
    return diag(*(Matrix(size,size,lambda i,j:
        eigenvalue if i==j else int(j==i+1)) for size in sizes))


def tensor_from_matrices(a0,a1,b0,b1):
    result={}
    for a,b,c,d in product(range(4),repeat=4):
        value=F(a0[a,c]*b0[b,d]+a1[a,c]*b1[b,d])
        if value:
            result[64*a+16*b+4*c+d]=value
    return result


def main():
    partitions=[(2,1,1),(2,2),(3,1),(4,)]
    table=[]
    for left in partitions:
        row=[]
        for right in partitions:
            a=jordan(left,1);b=jordan(right,-1)
            kernel=16-(eye(16)+kronecker_product(a,b)).rank()
            assert kernel==sum(min(p,q) for p in left for q in right)
            row.append(kernel)
        table.append(row)
    assert table==[[10,8,7,4],[8,8,6,4],[7,6,6,4],[4,4,4,4]]
    print('PASS all non-scalar single-eigenvalue Jordan partition pairs:',table)
    p=diag(1,0,0,0);q=eye(4)-p;n=zeros(4);n[0,1]=1
    cases={'semisimple':tensor_from_matrices(p,q,q,p),
           'nilpotent':tensor_from_matrices(n,eye(4),eye(4),n)}
    for name,w in cases.items():
        assert [flatten(w,(i,)) for i in range(4)]==[4]*4
        assert [flatten(w,cut) for cut in ((0,1),(1,2),(0,2))]==[6,6,2]
        dual=orthogonal_windows(reconstructed_windows(w))
        dimensions=[len(nullspace_sparse(insertion_columns(dual,i,3),64))
                    for i in range(4)]
        assert dimensions==[12]*4
        if name == 'nilpotent':
            words = [(0,0,1),(0,0,2),(0,0,3),(0,2,1),(0,3,1)]
            words += [(a,0,j) for a in (2,3) for j in (1,2,3)]
            basis = [{16*a+4*b+c:F(1)} for a,b,c in words]
            basis.append({0:F(1),5:F(1),17:F(1)})
            assert len(basis)==12 and exact_rank(basis)==12
            equations=list(insertion_columns(dual,0,3))
            assert all(sum(value*vector.get(i,F(0)) for i,value in row.items())==0
                       for row in equations for vector in basis)
        print(f'PASS {name}: singleton ranks 4, pair ranks (6,6,2), W3={dimensions}')
    print('PASS exact finite checks; full pencil coverage is proved in manuscript')


if __name__=='__main__':
    main()
