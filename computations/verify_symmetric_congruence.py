#!/usr/bin/env python3
"""Exact centroid/root checks for the symmetric tensor congruence lemma."""
from fractions import Fraction as F
from itertools import product
from sympy import Matrix, eye
from frobenius_tensor import digits, flatten, nullspace_sparse
from tensor_symmetrization import invertible_witness


def index(word):
    result=0
    for x in word:
        result=4*result+x
    return result


def act(tensor, matrix, factor):
    result={}
    for encoded,value in tensor.items():
        word=list(digits(encoded));source=word[factor]
        for target in range(4):
            coefficient=F(matrix[target,source])
            if coefficient:
                word[factor]=target
                key=index(word)
                result[key]=result.get(key,F(0))+coefficient*value
    return {i:x for i,x in result.items() if x}


def all_act(tensor, matrices):
    result=tensor
    for factor,matrix in enumerate(matrices):
        result=act(result,matrix,factor)
    return result


def centroid_dimension(tensor):
    equations={}
    for column in range(16):
        elementary=Matrix.zeros(4)
        elementary[column//4,column%4]=1
        base=act(tensor,elementary,0)
        for factor in range(1,4):
            other=act(tensor,elementary,factor)
            for key in base.keys() | other.keys():
                value=base.get(key,F(0))-other.get(key,F(0))
                if value:
                    equations.setdefault((factor,key),{})[column]=value
    return len(nullspace_sparse(equations.values(),16))


def main():
    # Semisimple centroid: a direct sum of four one-dimensional tensors.
    diagonal={index((i,)*4):F(1) for i in range(4)}
    assert [flatten(diagonal,(i,)) for i in range(4)]==[4]*4
    assert centroid_dimension(diagonal)==4
    matrices=[Matrix.diag(*((i+j+1)**4 for i in range(4))) for j in range(3)]
    root=Matrix.diag(*((i+1)*(i+2)*(i+3) for i in range(4)))
    assert all_act(diagonal,[eye(4)]+matrices)==all_act(diagonal,[root]*4)
    print('PASS semisimple centroid: dimension four, common fourth root')
    # Nonsemisimple centroid: trace of the local algebra k[t]/(t^4).
    local={index(word):F(1) for word in product(range(4),repeat=4) if sum(word)==3}
    assert [flatten(local,(i,)) for i in range(4)]==[4]*4
    assert centroid_dimension(local)==4
    shift=Matrix(4,4,lambda i,j:int(j==i+1))
    matrices=[eye(4)+a*shift for a in (1,2,3)]
    for matrix in matrices:
        assert all(act(local,matrix,i)==act(local,matrix,0) for i in range(4))
    total=matrices[0]*matrices[1]*matrices[2]
    nilpotent=total-eye(4)
    assert nilpotent**4==Matrix.zeros(4)
    root=eye(4)+nilpotent/4-3*nilpotent**2/32+7*nilpotent**3/128
    assert root**4==total
    assert all_act(local,[eye(4)]+matrices)==all_act(local,[root]*4)
    print('PASS nilpotent centroid: dimension four, exact polynomial fourth root')
    # Common changes of basis preserve the congruence equality.
    change=eye(4)
    change[0,1]=2;change[2,3]=-3
    lhs=all_act(local,[change]+[change*m for m in matrices])
    rhs=all_act(local,[change*root]*4)
    assert lhs==rhs
    print('PASS independent equivalence becomes one common congruence')
    basis=[Matrix(4,4,lambda a,b:int(a==b==i)) for i in range(4)]
    assert all(matrix.det()==0 for matrix in basis)
    assert invertible_witness(basis).det()!=0
    assert invertible_witness([basis[0],basis[1]]) is None
    assert nullspace_sparse([{0:F(0)}],2)==[{0:F(1)},{1:F(1)}]
    print('PASS determinant-span witnesses and all-singular rejection, zero coefficients')
    print('PASS exact representative checks; general lemma is proved in the manuscript')


if __name__=='__main__':
    main()
