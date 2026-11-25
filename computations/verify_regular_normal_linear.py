"""Rational checks of the marked PBW criterion and the unmarked elliptic charts."""

from itertools import product

import sympy as sp

from verify_quantum_triple_face import outer

TRIPLES = ((1, 2, 3), (2, 3, 1), (3, 1, 2))


def add(target, source, scale=1):
    for word, value in source.items():
        target[word] = target.get(word, 0) + scale * value


def verify():
    p, q, s = 1, 2, 4
    assert (p**3 + q**3 + s**3)**3 != 27*(p*q*s)**3
    relations = [{(b, c): p, (c, b): q, (a, a): s} for a, b, c in TRIPLES]
    coeff = sp.Matrix(3, 3, sp.symbols('m:9'))
    constants = sp.Matrix(sp.symbols('b:3'))
    overlap = {}
    defect = {}
    for a, _, _ in TRIPLES:
        axis = {(a,): 1}
        alpha = {(b,): coeff[a-1, b-1] for b in range(1, 4)}
        add(overlap, outer(axis, relations[a-1]))
        add(overlap, outer(relations[a-1], axis), -1)
        add(defect, outer(alpha, axis))
        add(defect, outer(axis, alpha), -1)
    assert all(v == 0 for v in overlap.values())
    assert all(defect[i, j] == coeff[j-1, i-1] - coeff[i-1, j-1]
               for i, j in product(range(1, 4), repeat=2))
    alpha = (coeff + coeff.T)/2
    rel_matrices = [sp.Matrix(3, 3, lambda i, j: r.get((i+1, j+1), 0)) for r in relations]

    def translate(linear, scalar, shift):
        correction = sp.Matrix.vstack(*[shift.T*(r+r.T) for r in rel_matrices])
        squares = sp.Matrix([(shift.T*r*shift)[0] for r in rel_matrices])
        return linear-correction, scalar+linear*shift-squares

    e = sp.Matrix(sp.symbols('e:3'))
    f = sp.Matrix(sp.symbols('f:3'))
    twice = translate(*translate(alpha, constants, e), f)
    once = translate(alpha, constants, e+f)
    assert all(sp.expand(x) == 0 for a, b in zip(twice, once) for x in a-b)
    normalized, _ = translate(alpha, constants, sp.Matrix([alpha[i, i]/(2*s) for i in range(3)]))
    assert normalized.is_symmetric() and all(normalized[i, i] == 0 for i in range(3))

    # The four-tensor is cyclic alternating for every symmetric linear correction.
    tensor = {}
    for a, _, _ in TRIPLES:
        lifted = dict(relations[a-1])
        for i in range(1, 4):
            for word in ((0, i), (i, 0)):
                lifted[word] = -alpha[a-1, i-1]/2
        lifted[0, 0] = -constants[a-1]
        central = {(0, a): 1, (a, 0): -1}
        add(tensor, outer(central, lifted))
        add(tensor, outer(lifted, central))
    assert all(sp.expand(tensor.get(w[1:]+w[:1], 0) + z) == 0 for w, z in tensor.items())


if __name__ == '__main__':
    verify()
    print('PASS: overlap, PBW symmetry, affine translation action, unique centered representative and cyclic tensor')
