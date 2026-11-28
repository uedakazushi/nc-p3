"""Certify the binary transfer, smooth secant sextic and unique triple octic."""
from itertools import combinations
from functools import reduce
from math import gcd, lcm

import numpy as np
from flint import nmod_mat, nmod_poly, fmpq_mat, fmpz_mat
from sympy import Matrix

from verify_quadric_biform_window import (
    PRIME, basis, coefficients, flip, forms, kernel, multiply, pivots, window,
)
from verify_quadric_kernel_projection import PAIRS, symmetric_root
from verify_opposite_root_triple_octics import monomials, octic_evaluation, quotient_projection


def quotient(v, a, b):
    old = [multiply(h,k) for h in v for k in basis(a-2,b-2)]
    full = basis(a,b)
    c = coefficients(old+full,a,b)
    positions = pivots(c)
    new = [full[k-len(old)] for k in positions[len(old):]]
    ann = kernel(c[:,:len(old)].T).T
    return new, fmpq_mat((ann@c[:,positions[len(old):]]).tolist()).inv()*fmpq_mat(ann.tolist())


def quartic_basis(v, da, db):
    z,_ = quotient(v,2,2)
    _,q = quotient(v,2+da,2+db)
    maps = [q*fmpq_mat(coefficients([multiply(h,k) for h in z],2+da,2+db).tolist())
            for k in basis(da,db)]
    denominator = lcm(*(int(e.denominator) for k in maps for e in k.entries()))
    first,second = [np.array([[int(k[j,i]*denominator) for j in range(4)] for i in range(5)],
                             dtype=object) for k in maps]
    eq = []
    for unit in np.eye(16,dtype=object):
        u = unit.reshape(4,4)
        eq.append((second@u[:,:3]-first@u[:,1:]).flatten())
    sol = kernel(np.array(eq,dtype=object).T)
    assert sol.shape == (16,1)
    u = sol.reshape(4,4)
    coords = np.column_stack((first@u,second@u[:,3]))
    return nmod_mat((coords%PRIME).tolist(),PRIME)


def sextic(v):
    p = PRIME
    transfer = quartic_basis(v,1,0).inv()*quartic_basis(v,0,1)
    powers = [(i,j) for i in range(7) for j in range(7-i)]
    rows,values = [],[]
    ruling = []
    for z in range(9):
        evaluation = nmod_mat([[z**k for k in range(5)]],p)*transfer
        vector = evaluation.tolist()[0]
        ruling.append([vector[:4],vector[1:]])
    for a,b in combinations(range(9),2):
        det = int(nmod_mat(ruling[a]+ruling[b],p).det())
        rows.append([pow(a+b,i,p)*pow(a*b,j,p)%p for i,j in powers])
        values.append([det*pow((a-b)**2,-1,p)%p])
    matrix = nmod_mat(rows,p)
    reduced,rank = matrix.transpose().rref()
    assert rank == len(powers)
    pos = [next(j for j in range(reduced.ncols()) if reduced[i,j]) for i in range(rank)]
    coeff = nmod_mat([rows[i] for i in pos],p).solve(nmod_mat([values[i] for i in pos],p))
    assert matrix*coeff == nmod_mat(values,p)
    return {(i,j,6-i-j):int(coeff[k,0]) for k,(i,j) in enumerate(powers)}


def smooth_and_diagonal(form):
    p = PRIME
    rows = []
    words = [(i,j,13-i-j) for i in range(14) for j in range(14-i)]
    for axis in range(3):
        for i in range(9):
            for j in range(9-i):
                mon = (i,j,8-i-j)
                row = []
                for word in words:
                    source = tuple(x-y+int(h==axis) for h,(x,y) in enumerate(zip(word,mon)))
                    row.append(source[axis]*form.get(source,0))
                rows.append(row)
    assert nmod_mat(rows,p).rank() == len(words)
    diag = [0]*13
    for (i,j,_),c in form.items():
        diag[i+2*j] += 2**i*c
    diagonal = nmod_poly(diag,p)
    assert diagonal.degree() == 12 and diagonal.gcd(diagonal.derivative()).degree() == 0


def verify_triple_ideal(edge):
    changes = symmetric_root()
    changes *= lcm(*(int(entry.q) for entry in changes))
    changes /= reduce(gcd,map(int,changes))
    rows = []
    slices = edge.reshape(4,4,6)
    for i in range(5):
        for j in range(5-i):
            for k in range(5-i-j):
                x = [1,i,j,k]
                cut = sum(a*z for a,z in zip(x,slices)).T
                row = [int(fmpz_mat(cut[list(h),:].tolist()).det()) for h in combinations(range(6),4)]
                coeff = (Matrix([[x[a]*x[b] for a,b in PAIRS]])*changes).tolist()[0]
                root = Matrix(4,4,lambda a,b: coeff[PAIRS.index(tuple(sorted((a,b))))])
                row.extend(int(root.extract(h,l).det()) for h in combinations(range(4),2)
                           for l in combinations(range(4),2))
                rows.append(row)
    # The 35 evaluations are unisolvent for quaternary quartics.
    assert fmpz_mat(rows).rank() == fmpz_mat([row[15:] for row in rows]).rank() == 15


def main():
    smooth_and_diagonal(sextic(forms()))
    edges,_,_ = window()
    verify_triple_ideal(edges[0])
    for r in (edges[1],flip(edges[2])):
        left = np.kron(r,np.eye(4,dtype=object))
        right = np.kron(np.eye(4,dtype=object),flip(r))
        assert kernel(np.column_stack((left,right))).shape[1] == 0
    rel = Matrix((edges[1]%PRIME).tolist())
    assert nmod_mat(octic_evaluation(rel).tolist(),PRIME).rank() == 164
    words,proj = quotient_projection(rel,4)
    multipliers = [nmod_mat(proj[:,[words[tuple(sorted((i,)+w))] for w in monomials(3)]].tolist(),PRIME)
                   for i in range(4)]
    affine = multipliers[0].inv()*multipliers[1]
    poly = affine.charpoly()
    assert poly.degree() == 20 and poly.gcd(poly.derivative()).degree() == 0
    print("PASS: smooth secant sextic, twelve simple diagonal intersections, twenty reduced marks, unique triple octic and both repeat obstructions")


if __name__ == "__main__":
    main()
