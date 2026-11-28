"""Binary cubic cone middle, its secant sextic, marked divisors and fold."""
from itertools import combinations, combinations_with_replacement, product
from math import lcm
from random import Random

import numpy as np
from flint import fmpq_mat, nmod_mat
from sympy import Matrix

from verify_opposite_root_triple_octics import monomials, octic_evaluation, quotient_projection
from verify_quadric_biform_window import PRIME, finite_projection, flip, kernel, pivots
from verify_quadric_kernel_projection import multiply
from verify_quadric_ruling_identifiability import smooth_and_diagonal
from verify_single_edge_instanton_extension import companion


def keys(degree):
    return [(i, j) for j in range(degree//2+1) for i in range(degree-2*j+1)]


def basis(degree):
    return [{key: 1} for key in keys(degree)]


def products(left, right):
    return [multiply(f, g) for f, g in product(left, right)]


def coefficients(forms, degree):
    return np.array([[f.get(key, 0) for f in forms] for key in keys(degree)], dtype=object)


def cone_edges():
    rng = Random(2867)
    v = [{key: rng.randint(-4, 4) for key in keys(4)} for _ in range(4)]
    assert len(pivots(coefficients(products(v, basis(2)), 6))) == 16
    columns = coefficients(products(v, basis(1)) + basis(5), 5)
    positions = pivots(columns)
    assert positions[:8] == list(range(8)) and len(positions) == 12
    new = [basis(5)[i-8] for i in positions[8:]]
    ann = kernel(columns[:, :8].T).T
    quotient = fmpq_mat((ann @ columns[:, positions[8:]]).tolist()).inv()*fmpq_mat(ann.tolist())
    den = lcm(*(int(t.denominator) for t in quotient.entries()))
    q = np.array([[int(quotient[i,j]*den) for j in range(12)] for i in range(4)], dtype=object)
    my = (q @ coefficients(products(basis(2), basis(3)), 5)).reshape(16, 6)
    square = [multiply(v[i], v[j]) for i,j in combinations_with_replacement(range(4), 2)]
    xm = kernel(coefficients(products(square, basis(1))+products(v, new), 9))[20:, :]
    assert xm.shape == my.shape == (16, 6)
    return xm, my


def field(matrix):
    return nmod_mat((matrix % PRIME).tolist(), PRIME)


def binary_coordinates(edge):
    first, second, third, weight_two = edge.reshape(4, 4, 6).transpose(1, 2, 0)
    equations = []
    for unit in np.eye(16, dtype=object):
        change = unit.reshape(4, 4)
        equations.append((second @ change[:, :3] - first @ change[:, 1:]).flatten())
    solutions = kernel(np.array(equations, dtype=object).T)
    assert solutions.shape == (16, 1)
    change = solutions[:, 0].reshape(4, 4)
    coords = np.column_stack((first @ change, second @ change[:, -1], third @ change[:, -1]))
    assert np.array_equal(coords[:, 1:5], second @ change)
    assert np.array_equal(coords[:, 2:6], third @ change)
    return field(change), field(coords).inv()*field(weight_two @ change)


def secant_sextic(transfer):
    powers = [(i,j) for i in range(7) for j in range(7-i)]
    ruling = []
    for z in range(9):
        quintic = nmod_mat([[z**i for i in range(6)]], PRIME)*transfer
        ruling.append([[z**i for i in range(4)], quintic.tolist()[0]])
    rows, values = [], []
    for a,b in combinations(range(9), 2):
        det = int(nmod_mat(ruling[a]+ruling[b], PRIME).det())
        rows.append([pow(a+b,i,PRIME)*pow(a*b,j,PRIME) % PRIME for i,j in powers])
        values.append([det*pow((a-b)**2, -1, PRIME) % PRIME])
    matrix = nmod_mat(rows, PRIME)
    red, rank = matrix.transpose().rref()
    assert rank == len(powers)
    pos = [next(j for j in range(red.ncols()) if red[i,j]) for i in range(rank)]
    coeff = nmod_mat([rows[i] for i in pos], PRIME).solve(nmod_mat([values[i] for i in pos], PRIME))
    assert matrix*coeff == nmod_mat(values, PRIME)
    return {(i,j,6-i-j): int(coeff[k,0]) for k,(i,j) in enumerate(powers)}


def squarefree_cubics(rel, change):
    words, proj = quotient_projection(rel, 4)
    maps = [nmod_mat(proj[:, [words[tuple(sorted((i,)+w))] for w in monomials(3)]].tolist(),
                     PRIME) for i in range(4)]
    ops = [maps[0].inv()*f for f in maps]
    poly = ops[1].charpoly()
    assert poly.degree() == 20 and poly.gcd(poly.derivative()).degree() == 0
    inv = change.inv()
    a,b,c,d = [sum((int(inv[i,j])*ops[j] for j in range(4)), 0*ops[0]) for i in range(4)]
    disc = b**2*c**2-4*a*c**3-4*b**3*d-27*a**2*d**2+18*a*b*c*d
    assert disc.rank() == 20


def fixed_tangent():
    cartan = np.zeros((12, 6), dtype=object)
    for i,j in product(range(4), range(3)):
        cartan[3*i+j, i+j] = 1
    b = kernel(cartan.T)
    pairs = list(combinations(range(6), 2))
    equations = np.array([
        [sum(b[3*u+i,k]*b[3*v+j,l]-b[3*v+j,k]*b[3*u+i,l] for u,v in ((a,d),(d,a)))
         for k,l in pairs]
        for a,d in combinations_with_replacement(range(4),2) for i,j in combinations(range(3),2)
    ], dtype=object)
    sol = kernel(equations)
    assert sol.shape == (15, 1)
    form = np.zeros((6, 6), dtype=object)
    for (i,j), val in zip(pairs, sol[:,0]):
        form[i,j], form[j,i] = val, -val
    assert len(pivots(form)) == 6


def main():
    fixed_tangent()
    xm, my = cone_edges()
    for edge in (xm, flip(xm), my, flip(my)):
        finite_projection(edge)
    change, transfer = binary_coordinates(my)
    smooth_and_diagonal(secant_sextic(transfer))
    rel = Matrix((my % PRIME).tolist())
    squarefree_cubics(rel, change)
    assert nmod_mat(octic_evaluation(rel).tolist(), PRIME).rank() == 164
    companion(Matrix(my.tolist()))
    print("PASS: four finite projections, smooth sextic, 20 separable marked cubics, unique octic, rank-ten fold and universal rank-six tangent")


if __name__ == "__main__":
    main()
