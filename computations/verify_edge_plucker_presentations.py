"""Recover the linear rows of an edge and test its Grassmannian deformation at fixed input."""

from itertools import combinations
from random import Random

import numpy as np
from flint import nmod_mat

from verify_configuration_intersections import exchange, kernel, modular_points, partner, rational_example


PAIRS = list(combinations(range(6), 2))
QUADS = list(combinations(range(6), 4))


def sign(sequence):
    return (-1)**sum(a>b for a,b in combinations(sequence,2))


def check(points, partners, prime):
    values = nmod_mat([[int(x)*int(y) % prime for x in p for y in q] for p,q in zip(points,partners)], prime)
    relations = kernel(values, prime)
    assert relations.ncols() == 6
    coefficients = np.array([[int(x) for x in row] for row in relations.tolist()]).reshape(4,4,6)
    rng = Random(591)
    evaluation, tangent, contractions, quadratic = [], [], [], []
    for _ in range(100):
        x = np.array([1, *(rng.randrange(prime) for _ in range(3))])
        edge = np.tensordot(x, coefficients, axes=(0,0)) % prime
        p = []
        for pair in PAIRS:
            others = tuple(i for i in range(6) if i not in pair)
            scalar = int(nmod_mat(edge[:,others].tolist(),prime).det())
            p.append(sign((*pair,*others))*scalar % prime)
        evaluation.append(p)
        quadratic.append([a*b % prime for i,a in enumerate(p) for b in p[i:]])
        skew = np.zeros((6,6),dtype=np.int64)
        for (i,j), v in zip(PAIRS,p):
            skew[i,j],skew[j,i] = v,-v
        assert not np.any((edge@skew)%prime)
        for row in skew:
            contractions.append(np.outer(row,x).flatten().tolist())
        for quad in QUADS:
            block = np.zeros((15,15),dtype=np.int64)
            for i,pair in enumerate(PAIRS):
                if set(pair) <= set(quad):
                    other = tuple(k for k in quad if k not in pair)
                    block[i,:] = sign((*pair,*other))*p[PAIRS.index(other)]*np.array(p) % prime
            tangent.append(block.flatten().tolist())
    assert nmod_mat(evaluation,prime).rank() == 15
    assert nmod_mat(contractions,prime).rank() == 20
    assert nmod_mat(quadratic,prime).rank() == 85
    assert nmod_mat(tangent,prime).rank() == 189


if __name__ == '__main__':
    prime = 2657
    tensor, _, _ = rational_example()
    a = modular_points(tensor, prime)
    b = partner(a, tensor, False, prime)
    c = partner(b, exchange(tensor), True, prime)
    d = partner(c, tensor.T, False, prime)
    for p,q in ((a,b), (a,d), (c,b), (c,d)):
        check(p,q,prime)
    print("All four presentations recover four linear rows, 85 quadratic products and tangent kernel dimension 36.")
