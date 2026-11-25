"""Exact transverse tangent bounds and simultaneous lifts over the whole quadric curve."""

from collections import defaultdict
from itertools import product

import sympy as sp

from verify_smooth_quadric_graphs import WORDS, volume

AXES = range(4)
PIVOT = [1, 2, 3, 4, 6, 7]
EXTRA = [j for j in range(16) if j not in PIVOT]
LETTER = ((1, 0), (-1, 0), (0, 1), (0, -1))
m = sp.symbols('m')


def charge(word):
    return tuple(sum(LETTER[j][k] for j in word) for k in range(2))


def index(word, side):
    a, b, c, d = word[side:]+word[:side]
    return 4*a+b, 4*c+d


def flatten(tensor, side):
    matrix = sp.zeros(16)
    for w, val in tensor.items():
        matrix[index(w, side)] = val
    return matrix


def tangent_letter(b):
    c, d = (2, 3) if b < 2 else (0, 1)
    return {(b,c,b,d):1, (b,d,b,c):-1, (c,b,d,b):-1, (d,b,c,b):1}


def unit_laurent(value):
    num, den = sp.factor(value).as_numer_denom()
    return value != 0 and all(poly in (m, m-1, m+1)
                              for part in (num, den) for poly, _ in sp.factor_list(part, m)[1])


def verify():
    base = volume(m)
    data = []
    eye = sp.eye(16)
    for side in range(2):
        mat = flatten(base, side)
        pivot = mat[PIVOT, PIVOT]
        assert sp.factor(pivot.det()) == -4*m**2*(m-1)**2*(m+1)**2
        inv = pivot.inv().applyfunc(sp.factor)
        left = (eye[EXTRA, :] - mat[EXTRA, PIVOT]*inv*eye[PIVOT, :]).applyfunc(sp.factor)
        right = (eye[:, EXTRA] - eye[:, PIVOT]*inv*mat[PIVOT, EXTRA]).applyfunc(sp.factor)
        data.append((inv, left, right))
    blocks = defaultdict(list)
    for word in WORDS:
        blocks[charge(word)].append(word)
    tangent = [tangent_letter(b) for b in AXES]
    for ch, words in blocks.items():
        rows = []
        for side, (_, left, right) in enumerate(data):
            for i, j in product(range(10), repeat=2):
                if charge(divmod(EXTRA[i], 4)+divmod(EXTRA[j], 4)) != ch:
                    continue
                rows.append([left[i, index(w, side)[0]]*right[index(w, side)[1], j] for w in words])
        linear = sp.Matrix(len(rows), len(words), sum(rows, []))
        sample = linear.subs(m, 2)
        cols = list(sample.rref()[1])
        pivots = list(sample[:, cols].T.rref()[1])
        assert unit_laurent(sp.factor(linear[pivots, cols].det(method='domain-ge')))
        orbit = []
        for pos, a, b in product(AXES, repeat=3):
            col = [base[w[:pos]+(b,)+w[pos+1:]] if w[pos] == a else 0 for w in words]
            if any(col):
                orbit.append(col)
        orbit = sp.Matrix(len(orbit), len(words), sum(orbit, [])).T
        cols2 = list(orbit.subs(m, 2).rref()[1])
        rows2 = list(orbit[:, cols2].subs(m, 2).T.rref()[1])
        assert unit_laurent(sp.factor(orbit[rows2, cols2].det(method='domain-ge')))
        bound = len(words)-len(cols)-len(cols2)
        axes = [i for i in AXES if ch == tuple(2*v for v in LETTER[i])]
        assert bound == (2 if ch == (0, 0) else len(axes))
        for i in axes:
            vec = sp.Matrix([tangent[i].get(w, 0) for w in words])
            assert all(sp.factor(v) == 0 for v in linear*vec)
            extra = orbit[:, cols2].row_join(vec)
            rr = list(extra.subs(m, 2).T.rref()[1])
            assert unit_laurent(sp.factor(extra[rr, :].det(method='domain-ge')))
    for side, (inv, left, right) in enumerate(data):
        mats = [flatten(t, side) for t in tangent]
        for a, b in [(i, i) for i in AXES] + list(product((0, 1), (2, 3))):
            value = left*mats[a][:, PIVOT]*inv*mats[b][PIVOT, :]*right
            value += left*mats[b][:, PIVOT]*inv*mats[a][PIVOT, :]*right
            assert all(sp.factor(v) == 0 for v in value)


if __name__ == '__main__':
    verify()
    print('PASS: all parameter-uniform tangent bounds, four transverse lines and the four exact attraction planes')
