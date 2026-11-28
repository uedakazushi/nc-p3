"""Recover the five-by-five symmetric linear block from the weighted quadric sample."""
from itertools import combinations, combinations_with_replacement, product
from random import Random

import numpy as np
from flint import nmod_mat
from sympy import Symbol, expand

from verify_quadric_biform_window import kernel, pivots
from verify_quadric_cone_middle import PRIME, basis, coefficients, cone_edges, field, keys, products
from verify_quadric_kernel_projection import multiply


def linear_block():
    rng = Random(2867)
    v = [{key: rng.randint(-4,4) for key in keys(4)} for _ in range(4)]
    pos = pivots(coefficients(v+basis(4),4))
    extras = [basis(4)[k-4] for k in pos[4:]]
    pairs = list(combinations_with_replacement(range(4),2))
    square = [multiply(v[i],v[j]) for i,j in pairs]
    relations = kernel(coefficients(square+products(v,extras),8))
    assert relations.shape == (30,5)
    blocks = relations[10:,:].reshape(4,5,5)
    eq = [[b[i,k]*int(l == j)-b[j,k]*int(l == i) for k,l in product(range(5),repeat=2)]
          for b in blocks for i,j in combinations(range(5),2)]
    solutions = kernel(np.array(eq,dtype=object))
    assert solutions.shape == (25,1)
    change = solutions[:,0].reshape(5,5)
    assert field(change).rank() == 5
    return [field(b @ change) for b in blocks]


def minors(matrix, pairs):
    return [int(nmod_mat([[matrix[i,j] for j in cols] for i in rows], PRIME).det())
            for rows,cols in pairs]


def main():
    blocks = linear_block()
    xm,_ = cone_edges()
    slices = [field(r.T) for r in xm.reshape(4,4,6)]
    sym = list(combinations(range(5),4))
    sym_pairs = [(sym[i],sym[j]) for i,j in combinations_with_replacement(range(5),2)]
    edge_pairs = [(rows,range(4)) for rows in combinations(range(6),4)]
    rows = []
    for i in range(5):
        for j in range(5-i):
            for k in range(5-i-j):
                x = [1,i,j,k]
                web = sum((u*b for u,b in zip(x,blocks)),0*blocks[0])
                edge = sum((u*b for u,b in zip(x,slices)),0*slices[0])
                rows.append(minors(web,sym_pairs)+minors(edge,edge_pairs))
    assert nmod_mat(rows,PRIME).rank() == 15
    assert nmod_mat([r[:15] for r in rows],PRIME).rank() == 15
    assert nmod_mat([r[15:] for r in rows],PRIME).rank() == 15
    t = Symbol("t")
    term = lambda n: (n+1)*(n+2)*(n+3)/6
    assert expand(term(t)-term(t-5)-5*term(t-6)+5*term(t-7)) == 20*t-45
    print("PASS: symmetric linear quintic block, matching fifteen marked quartics at the sample prime, genus-46 conductor")


if __name__ == "__main__":
    main()
