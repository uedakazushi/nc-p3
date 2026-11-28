"""Exact biform window, four finite edges, and an unbalanced simple representation."""
from functools import reduce
from itertools import combinations, combinations_with_replacement, product
from math import gcd, lcm

import numpy as np
from flint import fmpq_mat, fmpz_mat, nmod_mat

from verify_quadric_kernel_projection import basis, forms, monomials, multiply

PRIME = 1009


def coefficients(polys, a, b):
    return np.array([[poly.get(key,0) for poly in polys] for key in monomials(a,b)], dtype=object)


def primitive(matrix):
    for j in range(matrix.shape[1]):
        matrix[:,j] //= reduce(gcd, matrix[:,j].tolist())
    return matrix


def kernel(matrix):
    solution, count = fmpz_mat(matrix.tolist()).nullspace()
    columns = [[int(solution[i,j]) for j in range(count)] for i in range(solution.nrows())]
    return primitive(np.array(columns,dtype=object).reshape(solution.nrows(),count))


def pivots(matrix):
    result, rank = fmpq_mat(matrix.tolist()).rref()
    return [next(j for j in range(result.ncols()) if result[i,j]) for i in range(rank)]


def half_window(v, a, b):
    lower, full = basis(a-2,b-2), basis(a,b)
    old = [multiply(f,h) for f,h in product(v,lower)]
    all_columns = coefficients(old+full,a,b)
    positions = pivots(all_columns)
    assert positions[:8] == list(range(8)) and len(positions) == 12
    new = [full[k-8] for k in positions[8:]]
    annihilator = kernel(all_columns[:,:8].T).T
    transition = fmpq_mat((annihilator@all_columns[:,positions[8:]]).tolist())
    quotient = transition.inv()*fmpq_mat(annihilator.tolist())
    denominator = lcm(*(int(t.denominator) for t in quotient.entries()))
    quotient = np.array([[int(denominator*quotient[i,j]) for j in range(12)] for i in range(4)],dtype=object)
    products = [multiply(y,h) for y,h in product(basis(1,1),basis(a-1,b-1))]
    middle_y = (quotient@coefficients(products,a,b)).reshape(16,6)
    old_square = [multiply(multiply(v[i],v[j]),h) for i,j in combinations_with_replacement(range(4),2)
                  for h in lower]
    products = [multiply(f,h) for f,h in product(v,new)]
    full_relations = kernel(coefficients(old_square+products,a+2,b+2))
    x_middle = full_relations[20:,:]
    assert x_middle.shape == (16,6) and len(pivots(x_middle)) == 6
    da,db = (0,1) if a == 2 else (1,0)
    higher = [multiply(f,h) for f,h in product(v,basis(a+da-2,b+db-2))]
    dual = kernel(coefficients(higher,a+da,b+db).T).T
    assert dual.shape == (3,15)
    arrows = []
    for h in new:
        evaluation = dual@coefficients([multiply(h,l) for l in basis(da,db)],a+da,b+db)
        arrows.append(np.stack([evaluation[:,1],-evaluation[:,0]]))
    return x_middle, middle_y, arrows


def flip(edge):
    return edge.reshape(4,4,6).transpose(1,0,2).reshape(16,6)


def finite_projection(edge):
    matrices = np.array(edge%PRIME,dtype=np.int64).reshape(4,4,6)
    values = []
    for i,j in monomials(4,4):
        if i+j > 4:
            continue
        point = [1,i,j,2+3*i+5*j]
        cut = nmod_mat(sum(c*z for c,z in zip(point,matrices)).T.tolist(), PRIME)
        values.append([int(nmod_mat([[cut[r,h] for h in range(4)] for r in indices],PRIME).det())
                       for indices in combinations(range(6),4)])
    assert nmod_mat(values,PRIME).rank() == 15


def finite_frobenius(relations):
    frames = relations
    for target in (4,1,0):
        new = []
        for i in range(4):
            left = np.kron(frames[i],np.eye(4,dtype=object))
            right = np.kron(np.eye(4,dtype=object),frames[(i+1)%4])
            common = kernel(np.column_stack((left,right)))
            assert common.shape[1] == target
            new.append(primitive(left@common[:left.shape[1],:]))
        frames = new
        if target == 1:
            tensors = [v[:,0].reshape((4,)*4) for v in frames]
            top = tensors[0]
            for i,t in enumerate(tensors):
                cycle = top.transpose(tuple((i+j)%4 for j in range(4)))
                assert np.array_equal(t,cycle) or np.array_equal(t,-cycle)
                for n,r in ((1,4),(2,6),(3,4)):
                    assert len(pivots(t.reshape(4**n,-1))) == r
            assert len(pivots(top.transpose(0,2,1,3).reshape(16,16))) == 16


def maximal_ruling_ranks(v):
    for transpose in (False,True):
        rows = []
        for z in range(7):
            block = [[sum(h.get((j,i) if transpose else (i,j),0)*z**j for j in range(3))
                      for h in v] for i in range(3)]
            rows.append([int(fmpz_mat([[row[k] for k in indices] for row in block]).det())
                         for indices in combinations(range(4),3)])
        # Four binary sextics span a basepoint-free system iff suitable multiples fill some degree.
        powers = nmod_mat([[z**j for j in range(7)] for z in range(7)],PRIME)
        coeff = (powers.inv()*nmod_mat(rows,PRIME)).tolist()
        multiples = [[int(coeff[j-i][k]) if 0 <= j-i <= 6 else 0 for j in range(12)]
                     for k in range(4) for i in range(6)]
        assert nmod_mat(multiples,PRIME).rank() == 12


def compose(first, second, relation):
    return sum(relation[4*i+j]*(second[j]@first[i]) for i,j in product(range(4),repeat=2))


def simple_module(relations, middle, last):
    y = [np.eye(4,dtype=object)[k].reshape(2,2) for k in range(4)]
    for j in range(6):
        assert not np.any(compose(middle,y,relations[1][:,j]))
        assert not np.any(compose(y,last,relations[2][:,j]))
    equations = []
    for k in range(36):
        first = np.eye(36,dtype=object)[k].reshape(4,3,3)
        column = []
        for j in range(6):
            column.extend(compose(first,middle,relations[0][:,j]).flat)
            column.extend(compose(last,first,relations[3][:,j]).flat)
        equations.append(column)
    solutions = kernel(np.array(equations,dtype=object).T)
    assert solutions.shape == (36,4)
    # Choose a solution with zero common kernel and full total image.
    for solution in solutions.T:
        first = solution.reshape(4,3,3)
        if len(pivots(np.column_stack(first))) == len(pivots(np.vstack(first))) == 3:
            return
    raise AssertionError("No full corner action found")


def window():
    v = forms()
    first, second, middle = half_window(v,2,3)
    reverse_first, reverse_second, reverse_middle = half_window(v,3,2)
    return [first,second,flip(reverse_second),flip(reverse_first)], middle, [t.T for t in reverse_middle]


if __name__ == "__main__":
    edges, middle, last = window()
    for edge in edges:
        finite_projection(edge)
        finite_projection(flip(edge))
    finite_frobenius(edges)
    maximal_ruling_ranks(forms())
    simple_module(edges,middle,last)
    print("PASS: four finite edges, rational 14641 Frobenius tensor of opposite rank 16, simple of dimensions (3,3,2,2)")
