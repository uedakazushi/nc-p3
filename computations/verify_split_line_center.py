"""Quartic central sections of split line windows, including all root-of-unity jumps."""
from collections import Counter
from itertools import product

import sympy as sp

from flint import nmod_mat

from verify_isolated_split_component import PRIME, WORDS, residue, tangent_dimensions
from verify_split_spectral_line import A, I, line_tensor, multiply, torus

T = sp.Symbol('t')
LETTERS = range(4)
UNITS = (1, sp.I, -1, -sp.I)


def weight(value):
    value = sp.cancel(value)
    exponent = int(value.as_powers_dict().get(T, 0))
    return exponent, UNITS.index(sp.cancel(value/T**exponent))


def add(first, second):
    return first[0]+second[0], (first[1]+second[1])%4


def minus(w):
    return -w[0], -w[1]%4


def relations(w):
    pairs = tuple(product(LETTERS, repeat=2))
    rewrites = []
    for p in LETTERS:
        def value(q):
            return w.get(q[-p:]+q[:-p] if p else q, 0)
        matrix = sp.Matrix([[value(a+b) for a in pairs] for b in pairs]).rref()[0]
        moves = {}
        for j in range(6):
            (a, u), (b, v) = [(q, matrix[j,k]) for k, q in enumerate(pairs) if matrix[j,k]]
            factor = weight(-v/u)
            moves[a] = b, factor
            moves[b] = a, minus(factor)
        rewrites.append(moves)
    return rewrites


def paths(moves, p, length):
    # Store word = t^n i^m times its component's first word.
    images = {}
    reps = []
    for root in product(LETTERS, repeat=length):
        if root in images:
            continue
        reps.append(root)
        images[root] = root, (0, 0)
        queue = [root]
        for q in queue:
            for i in range(length-1):
                step = moves[(p+i)%4].get(q[i:i+2])
                if step is None:
                    continue
                pair, factor = step
                new = q[:i]+pair+q[i+2:]
                coeff = add(images[q][1], minus(factor))
                if new in images:
                    assert images[new] == (root, coeff)
                else:
                    images[new] = root, coeff
                    queue.append(new)
    assert len(reps) == sp.binomial(length+3, 3)
    return reps, images


def obstruction(graph, killed):
    seen = set()
    divisors = []
    for root in graph:
        if root in seen:
            continue
        weights = {root:(0, 0)}
        queue = [root]
        residuals = set()
        for v in queue:
            for u, factor in graph[v]:
                tentative = add(weights[v], factor)
                if u in weights:
                    residuals.add(add(tentative, minus(weights[u])))
                else:
                    weights[u] = tentative
                    queue.append(u)
        seen.update(weights)
        if killed.intersection(weights):
            continue
        gcd = sp.Poly(0, T, extension=sp.I)
        for n, m in residuals:
            polynomial = UNITS[m]*T**max(0,n)-T**max(0,-n)
            gcd = gcd.gcd(sp.Poly(polynomial, T, extension=sp.I))
        divisors.append(gcd.monic().as_expr())
    return Counter(divisors)


def quartic_center(w):
    moves = relations(w)
    reps, _ = zip(*(paths(moves, p, 4) for p in LETTERS))
    _, images = zip(*(paths(moves, p, 5) for p in LETTERS))
    graph = {(p,q):[] for p in LETTERS for q in reps[p]}
    killed = set()
    for p, a in product(LETTERS, repeat=2):
        nxt = (p+1)%4
        left = {images[p][q+(a,)][0]:((p,q), images[p][q+(a,)][1]) for q in reps[p]}
        right = {images[p][(a,)+q][0]:((nxt,q), images[p][(a,)+q][1]) for q in reps[nxt]}
        assert len(left) == len(right) == 35
        for key in left.keys() | right.keys():
            if key in left and key in right:
                (u, c), (v, d) = left[key], right[key]
                delta = add(c, minus(d))
                graph[u].append((v, delta))
                graph[v].append((u, minus(delta)))
            else:
                killed.add((left.get(key) or right[key])[0])
    return obstruction(graph, killed)


def verify_unipotent_axis(d):
    skew = sp.Matrix([[0, -1], [1, 0]])
    c = sp.Matrix([[1, 1], [0, 1]])*skew.inv()
    w = line_tensor(sp.Integer(1), (c, d, c.inv().T, d.inv().T))
    composite = I
    for i in LETTERS:
        def value(q):
            return w.get(q[-i:]+q[:-i] if i else q, 0)
        slices = [sp.Matrix(2, 2, lambda a,b:value((a,b,f,g)))
                  for f,g in product((2,3), repeat=2)]
        coeff = next(m for m in slices if m != sp.zeros(2))
        composite = skew*coeff.T*composite
    assert (composite/composite[0,0]).applyfunc(sp.cancel) == sp.Matrix([[1,0],[-2,1]])


def verify():
    singleton = line_tensor(T, (I,)*4)
    mixed = line_tensor(sp.I, multiply(torus(T, sp.Integer(1)), A))
    assert quartic_center(singleton) == Counter({T**2-1:6, 0:3, T**4-1:2})
    assert quartic_center(mixed) == Counter({1:14})
    skew = sp.Matrix([[0, -1], [1, 0]])
    for matrix in (sp.diag(T,1), sp.Matrix([[1,1],[0,1]])):
        verify_unipotent_axis(matrix*skew.inv())
    w = {q:residue(sp.sympify(v).subs(T,3)) for q,v in mixed.items()}
    assert tangent_dimensions(w) == [0,16,28,16,0]
    columns = list(product(LETTERS, repeat=3))
    orbit = nmod_mat([[w.get(q[:i]+(b,)+q[i+1:],0) if q[i]==a else 0
                       for i,a,b in columns] for q in WORDS], PRIME)
    assert orbit.rank() == 59
    print('All quartic-center jumps, the Jordan axis, tangent dimension 60 and orbit dimension 59 verified.')


if __name__ == '__main__':
    verify()
