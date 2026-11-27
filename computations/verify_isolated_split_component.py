"""Lift the exceptional split tensor exactly and bound its unrestricted tangent space."""
from itertools import product

import sympy as sp
from flint import nmod_mat

from verify_split_inverse_windows import BITS, PAIRS, bell, borders, positions

PRIME = 109
ROOT = 17
SQRT_MINUS_ONE = 33
WORDS = tuple(product(range(4), repeat=4))


def polynomials():
    d = sp.Symbol('d')
    return d, (d**2 + 1)**3 - 4*d**3


def reduce_polynomial(polynomial):
    d, modulus = polynomials()
    return sp.rem(polynomial, modulus, d, domain=sp.QQ_I)


def residue(value):
    d, _ = polynomials()
    real, imag = sp.expand(value.subs(d, ROOT) if hasattr(value, 'subs') else value).as_real_imag()
    a = sp.Rational(real + SQRT_MINUS_ONE*imag)
    return int(a.p) * pow(int(a.q), -1, PRIME) % PRIME


def projective_inverses(tensor):
    tensors = []
    for start in (0, 1):
        flat = sp.Matrix(4, 4, lambda i, j: tensor[positions(i,j,start)])
        # These cuts split into two binary blocks; adjugate entries are cubics.
        adj = flat.adjugate(method='berkowitz').T
        tensors.append({positions(i,j,start): reduce_polynomial(adj[i,j])
                        for i,j in product(range(4), repeat=2)})
    return tensors


def carry(tensor, maps, pairing):
    return {q: reduce_polynomial(sum(v*sp.prod(maps[i][q[pairing[i]],j[i]] for i in range(4))
                                     for j,v in tensor.items()))
            for q in BITS}


def spectrum_at(d):
    return [sp.Integer(1), (d**5-d**4+3*d**3-7*d**2+6*d-2)/2,
            (d**3+d**2+3*d-1)/2, d]


def exceptional_maps():
    return ([sp.eye(2)]*4, [sp.Matrix([[1,1],[-sp.I,sp.I]]),
                            sp.Matrix([[1,1],[sp.I,-sp.I]])/2,
                            sp.Matrix([[1,-1],[sp.I,sp.I]]),
                            sp.Matrix([[1,-1],[-sp.I,-sp.I]])/2])


def certificate():
    d, modulus = polynomials()
    tensor = bell(spectrum_at(d))
    maps = exceptional_maps()
    adjs = projective_inverses(tensor)
    transported = [carry(a,m,p) for a,m,p in zip(adjs,maps,PAIRS)]
    pivot = next(q for q in BITS if residue(transported[0][q]) and residue(transported[1][q]))
    for q in BITS:
        assert reduce_polynomial(transported[0][pivot]*transported[1][q]
                                 - transported[1][pivot]*transported[0][q]) == 0
    assert pow(SQRT_MINUS_ONE,2,PRIME) == PRIME-1
    assert modulus.subs(d,ROOT) % PRIME == 0
    assert residue(modulus.diff(d)) != 0
    # Unscaled inverses here give actual matrix Schur complements, not just
    # the projective equality verified over the characteristic-zero ring.
    duals = []
    for start in (0, 1):
        flat = nmod_mat([[residue(tensor[positions(i,j,start)]) for j in range(4)]
                         for i in range(4)], PRIME)
        inv = flat.inv().transpose()
        duals.append({positions(i,j,start): int(inv[i,j]) for i,j in product(range(4),repeat=2)})
    targets = [{q: residue(value) for q,value in carry(m,ls,pr).items()}
               for m,ls,pr in zip(duals,maps,PAIRS)]
    scale = targets[0][pivot] * pow(targets[1][pivot],-1,PRIME) % PRIME
    assert all((targets[0][q]-scale*targets[1][q]) % PRIME == 0 for q in BITS)
    slots = (0,1,0,1)
    def word(q, colors):
        return tuple(j+2*k for j,k in zip(q,colors))
    w = {word(q,slots):residue(v) for q,v in tensor.items()}
    w.update({word(q,[1-k for k in slots]):v for q,v in targets[0].items()})
    for i in range(2):
        first = {i,i+1}
        for q in BITS:
            values = borders(maps[i], PAIRS[i], i, q)
            for k in (0,1):
                colors = [slots[j] ^ (int(j not in first) ^ k) for j in range(4)]
                w[word(q,colors)] = residue(values[k]) * (scale if i == 1 and k == 0 else 1) % PRIME
    return w


def tangent_dimensions(w):
    equations = []
    for i in (0,1):
        flat = nmod_mat(16,16,[w.get(q[-i:]+q[:-i] if i else q,0) for q in WORDS],PRIME)
        left,n = flat.transpose().nullspace()
        right,m = flat.nullspace()
        assert n == m == 10
        equations.extend([[int(left[4*q[i]+q[(i+1)%4],j]*right[4*q[(i+2)%4]+q[(i+3)%4],k])
                           for q in WORDS] for j,k in product(range(n), range(m))])
    dims = []
    for n in range(5):
        cols = [j for j,q in enumerate(WORDS) if sum(t < 2 for t in q) == n]
        restricted = nmod_mat([[row[j] for j in cols] for row in equations], PRIME)
        dims.append(len(cols)-restricted.rank())
    assert nmod_mat(equations,PRIME).rank() == 196
    return dims


if __name__ == '__main__':
    dims = tangent_dimensions(certificate())
    assert dims == [0,16,28,16,0]
    print("Exact exceptional transporter over Q(i)[d]/((d^2+1)^3-4d^3).")
    print("Unrestricted affine tangent dimensions by normal weight:", dims)
