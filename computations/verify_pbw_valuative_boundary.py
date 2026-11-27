"""Determinants in descending Koszul bases and the full mixed-flip orbit rank."""
from itertools import combinations, product

import sympy as sp
from flint import nmod_mat

from verify_split_line_center import LETTERS, PRIME, T, WORDS, relations, residue, weight
from verify_split_spectral_line import A, line_tensor, multiply, torus


def labels(i, letter):
    return letter ^ (int(i%4 >= 2) if letter < 2 else i%2)


def desc(i, q):
    return all(labels(i+j, a) > labels(i+j+1, b) for j, (a,b) in enumerate(zip(q,q[1:])))


def pairing_orders(w):
    rewrites = relations(w)
    orders = []
    for i in LETTERS:
        descending = [q for q in rewrites[i] if desc(i,q)]
        assert len(descending) == 6
        for q in descending:
            target = rewrites[i][q][0]
            before = [min(2,labels(i+j,v)) for j,v in enumerate(q)]
            after = [min(2,labels(i+j,v)) for j,v in enumerate(target)]
            assert after == before[::-1]
            assert labels(i,target[0]) < labels(i+1,target[1])
        def coeff(q):
            return w.get(q[-i:]+q[:-i] if i else q,0)
        pivot = next(coeff(q) for q in WORDS if desc(i,q))
        triples = [q for q in product(LETTERS,repeat=3) if desc(i+1,q)]
        pair = sp.Matrix([[coeff((a,)+q)/pivot for q in triples] for a in LETTERS])
        orders.append(weight(pair.det())[0])
    assert orders == [2,1,0,-1]
    return orders


def pivots(matrix):
    reduced, rank = matrix.rref()
    return [next(j for j in range(matrix.ncols()) if reduced[i,j]) for i in range(rank)]


def orbit_rank(w):
    columns = tuple(product(LETTERS,repeat=3))
    jacobian = sp.Matrix([[w.get(q[:i]+(b,)+q[i+1:],0) if q[i]==a else 0
                           for i,a,b in columns] for q in WORDS])
    def sample(at):
        return nmod_mat([[residue(sp.sympify(v).subs(T,at)) for v in row]
                         for row in jacobian.tolist()], PRIME)
    matrix = sample(3)
    cols = pivots(matrix)
    rows = pivots(nmod_mat([[int(matrix[i,j]) for j in cols] for i in range(matrix.nrows())],
                           PRIME).transpose())
    assert len(cols) == len(rows) == 59
    det = jacobian.extract(rows,cols).det(method='domain-ge')
    assert sp.cancel(det+64*T**13*(T**2+1)) == 0
    for value in (sp.I, -sp.I):
        assert sample(value).rank() == 59


def quantum_pairing():
    # Coefficients are valuations, with one independent symbol for every interchange.
    weights = dict(zip(combinations(LETTERS,2), sp.symbols('c:6')))
    orders = [sum(weights[j,k] for k in LETTERS if k>j) for j in LETTERS]
    assert sum(orders) == sum(weights.values())


if __name__ == '__main__':
    tensor = line_tensor(sp.I, multiply(torus(T, sp.Integer(1)), A))
    print('Orders of first Koszul pairings', pairing_orders(tensor))
    orbit_rank(tensor)
    quantum_pairing()
    print('The affine orbit has rank 59 for every nonzero parameter.')
