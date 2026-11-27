"""Polynomial normal forms and stabilizer action in the exceptional mixed-flip fiber."""
from collections import defaultdict
from itertools import product

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from verify_isolated_split_component import WORDS
from verify_mixed_flip_slice import degree
from verify_split_spectral_line import A, line_tensor, multiply, torus

WEIGHTS = ((-1, 0), (-1, 2), (0, -1), (0, 2), (1, 1), (2, -1), (2, 0))
CHAMBERS = ((-3, -1), (-1, -3), (-1, 2), (2, -1), (1, 1))
UNITS = (1, sp.I, -1, -sp.I)


def pivots(matrix):
    return DomainMatrix.from_Matrix(matrix).to_field().rref()[1]


def cut(tensor, start):
    return sp.ImmutableMatrix(16, 16, [
        tensor.get(q[-start:] + q[:-start] if start else q, 0) for q in WORDS
    ])


def positive_support(weight):
    return tuple(h for h in WEIGHTS if sum(a*b for a, b in zip(h, weight)) > 0)


class TransverseSlice:
    def __init__(self):
        self.base = line_tensor(sp.I, multiply(torus(sp.I, sp.Integer(1)), A))
        self.words = defaultdict(list)
        self.rows = defaultdict(list)
        self.cuts = []
        self.basis = {}
        self.solver = {}
        for q in WORDS:
            self.words[degree(q)].append(q)
        for start in (0, 1):
            matrix = cut(self.base, start)
            left = sp.Matrix.hstack(*matrix.T.nullspace()).T
            right = sp.Matrix.hstack(*matrix.nullspace())
            rows, cols = pivots(matrix.T), pivots(matrix)
            inverse = matrix.extract(rows, cols).inv()
            g = sp.zeros(16)
            for j, k in product(range(6), repeat=2):
                g[cols[j], rows[k]] = inverse[j, k]
            assert matrix * g * matrix == matrix
            self.cuts.append((g, left, right))
            for j, k in product(range(10), repeat=2):
                row = {}
                for q in WORDS:
                    t = q[start:] + q[:start]
                    entry = left[j, 4*t[0]+t[1]] * right[4*t[2]+t[3], k]
                    if entry:
                        row[q] = entry
                h = degree(next(iter(row)))
                assert all(degree(q) == h for q in row)
                self.rows[h].append((start, j, k, row))
        for h, words in self.words.items():
            matrix = sp.Matrix([[row.get(q, 0) for q in words] for *_, row in self.rows[h]])
            ix = pivots(matrix.T)
            linear = matrix[list(ix), :]
            star = linear.conjugate().T
            self.solver[h] = (ix, star*(linear*star).inv()) if ix else ((), sp.zeros(len(words), 0))
            orbit = sp.Matrix([
                [self.base.get(q[:i]+(b,)+q[i+1:], 0) if q[i] == a else 0
                 for i, a, b in product(range(4), repeat=3)] for q in words
            ])
            transverse = matrix.col_join(orbit.conjugate().T).nullspace()
            if transverse:
                assert len(transverse) == 1
                self.basis[h] = dict(zip(words, transverse[0]))
        assert set(self.basis) == set(WEIGHTS) | {(0, 0)}

    def cone(self, weight):
        """Return monomial coefficients of the polynomial tensor, with constant term omitted."""
        def value(h):
            return sum(a*b for a, b in zip(weight, h))
        hs = positive_support(weight)
        maximum = max(map(value, self.words))
        ds = list(map(value, hs))
        exponents = [k for k in product(*(range(maximum//d+1) for d in ds))
                     if 0 < sum(a*d for a, d in zip(k, ds)) <= maximum]
        tensors, effective, products = {}, {}, {}
        for k in sorted(exponents, key=sum):
            h = tuple(sum(a*w[i] for a, w in zip(k, hs)) for i in (0, 1))
            if h not in self.words:
                continue
            rhs = [sp.zeros(16), sp.zeros(16)]
            for u, mats in products.items():
                v = tuple(a-b for a, b in zip(k, u))
                if v in effective:
                    for i in (0, 1):
                        rhs[i] += mats[i]*effective[v][i]
            if sum(k) == 1:
                tensors[k] = self.basis[h]
            else:
                small = [left*rhs[i]*right for i, (_, left, right) in enumerate(self.cuts)]
                out = [small[i][j, t] for i, j, t, _ in self.rows[h]]
                ix, solver = self.solver[h]
                vector = solver*sp.Matrix(len(ix), 1, [out[j] for j in ix])
                tensors[k] = {q: sp.expand(c) for q, c in zip(self.words[h], vector) if c}
            matrices = [cut(tensors[k], i) for i in (0, 1)]
            effective[k] = [matrices[i]-rhs[i] for i in (0, 1)]
            products[k] = [matrices[i]*self.cuts[i][0] for i in (0, 1)]
            for i, (_, left, right) in enumerate(self.cuts):
                assert left*effective[k][i]*right == sp.zeros(10)
        return hs, tensors


def actions(slice_):
    words = sorted(slice_.base)
    incidence = sp.Matrix([[int(q[i] == a) for i, a in product(range(4), repeat=2)]
                           for q in words])
    diagonal, left, right = [a.to_Matrix() for a in smith_normal_decomp(
        DomainMatrix.from_Matrix(incidence).convert_to(sp.ZZ))]
    rank = incidence.rank()
    assert list(diagonal.diagonal()[:rank]) == [1]*10 + [2]

    def target(q, rotation, flips):
        t = q[rotation:] + q[:rotation]
        return tuple(a ^ flips[2*i+a//2] for i, a in enumerate(t))

    def evaluate(rotation, flips, phases):
        result = []
        for h in WEIGHTS:
            img = {}
            for q, b in slice_.basis[h].items():
                if b:
                    power = sum(phases[4*i+j] for i, j in enumerate(q)) % 1
                    factor = sp.expand_complex(sp.exp(2*sp.pi*sp.I*power))
                    img[target(q, rotation, flips)] = sp.simplify(factor*b)
            q = next(iter(img))
            new = degree(q)
            basis = slice_.basis[new]
            factor = sp.simplify(img[q]/basis[q])
            assert all(sp.simplify(img.get(t, 0)-a*factor) == 0 for t, a in basis.items())
            result.append((WEIGHTS.index(new), factor))
        return result

    identity = list(range(7))
    reflected = [2, 5, 0, 6, 4, 1, 3]
    expected = [(identity, [1]*7), (reflected, [sp.I, sp.I, sp.I, -sp.I, 1, sp.I, sp.I])]
    for m in (0, 1):
        possibilities = []
        for flips in product(range(2), repeat=8):
            targets = [target(q, m, flips) for q in words]
            if set(targets) != set(words):
                continue
            y = left*sp.Matrix([UNITS.index(sp.cancel(slice_.base[t]/slice_.base[q]))
                                for q, t in zip(words, targets)])/4
            if any(not y[i].is_Integer for i in range(rank, len(y))):
                continue
            z = sp.zeros(16, 1)
            for i in range(rank):
                z[i] = y[i]/diagonal[i, i]
            possibilities.append(evaluate(m, flips, right*z))
            if m:
                break
        candidates = expected if m == 0 else [(identity, [1, -1, -sp.I, -1, -sp.I, -sp.I, 1])]
        assert possibilities == [list(zip(p, factors)) for p, factors in candidates]
    sign = evaluate(0, (0,)*8, right[:,rank-1]/2)
    assert sign == list(zip(identity, [-1, -1, -1, 1, 1, -1, 1]))


if __name__ == '__main__':
    slice_ = TransverseSlice()
    for chamber in CHAMBERS:
        slice_.cone(chamber)
    actions(slice_)
    print('All five polynomial attracting spaces, their stabilizer generators and the shift verified.')
