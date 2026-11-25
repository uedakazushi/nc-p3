"""Check the universal quantum chambers and exact Schur reconstruction over Q."""

from collections import Counter, defaultdict
from itertools import combinations, permutations, product

import sympy as sp

from verify_quantum_attraction_fiber import cross, dot
from verify_quantum_triple_face import basis_tensors, weight

AXES = tuple(range(4))
WORDS = tuple(product(AXES, repeat=4))
PAIRS = tuple(product(AXES, repeat=2))
PIVOT = tuple(i for i, (a, b) in enumerate(PAIRS) if a < b)
EXTRA = tuple(i for i in range(16) if i not in PIVOT)
CHAMBERS = ((0, 1, 4, 6), (0, 3, 5, 9), (0, 6, 8, 9))


def verify_chambers():
    chars = sorted({weight(w) for w in WORDS} - {(0,) * 4})
    walls = [w[:3] for w in chars]
    supports = set()
    for a, b in permutations(walls, 2):
        ray = cross(a, b)
        if not any(ray):
            continue
        edge = cross(a, ray)
        flags = [(dot(w, ray), dot(w, edge), dot(w, a)) for w in walls]
        for signs in product((-1, 1), repeat=3):
            supports.add(frozenset(i for i, flag in enumerate(flags)
                                   if next(s*v for s, v in zip(signs, flag) if v) > 0))
    maximal = {s for s in supports if not any(s < t for t in supports)}
    assert len(supports) == 336
    assert Counter(map(len, maximal)) == {17: 48, 20: 24}
    def support(cocharacter):
        return frozenset(i for i, w in enumerate(chars) if dot(cocharacter, w) > 0)
    assert maximal == {support(p) for c in CHAMBERS for p in permutations(c)}
    assert [max(dot(w, c) for w in chars) for c in CHAMBERS] == [13, 19, 13]


def cut(word, side):
    i, j, k, l = word[side:] + word[:side]
    return 4*i + j, 4*k + l


def flatten(tensor, side):
    matrix = sp.zeros(16)
    for word, value in tensor.items():
        matrix[cut(word, side)] = value
    return matrix


class CutEquations:
    def __init__(self, background, side):
        base = flatten(background, side)
        inv = base.extract(PIVOT, PIVOT).inv()
        self.base_determinant = base.extract(PIVOT, PIVOT).det()
        self.side = side
        select = sp.eye(16)
        self.left = select[EXTRA, :] - base.extract(EXTRA, PIVOT)*inv*select[PIVOT, :]
        self.right = select[:, EXTRA] - select[:, PIVOT]*inv*base.extract(PIVOT, EXTRA)

    def residual(self, tensor):
        matrix = flatten(tensor, self.side)
        core = matrix.extract(PIVOT, PIVOT)
        assert core.det() == self.base_determinant
        return matrix.extract(EXTRA, EXTRA) - matrix.extract(EXTRA, PIVOT)*core.inv()*matrix.extract(PIVOT, EXTRA)


def orbit_complement(background, words):
    # At the rational backgrounds tested the orthogonal complement is invariant
    # under the stabilizer (unitary form at q=1, common torus and cycles otherwise).
    orbit = []
    for side in AXES:
        for a, b in permutations(AXES, 2):
            column = [background.get(w[:side] + (b,) + w[side+1:], 0)
                      if w[side] == a else 0 for w in words]
            if any(column):
                orbit.append(column)
    if not orbit:
        return sp.eye(len(words))
    return sp.Matrix.hstack(*sp.Matrix(orbit).nullspace())


class WeightRecursion:
    def __init__(self, background, cocharacter):
        self.background = background
        self.cuts = [CutEquations(background, side) for side in range(2)]
        words_by_char = defaultdict(list)
        for w in WORDS:
            if dot(weight(w), cocharacter) > 0:
                words_by_char[weight(w)].append(w)
        self.blocks = []
        for char in sorted(words_by_char, key=lambda w: (dot(w, cocharacter), w)):
            words = words_by_char[char]
            normal = orbit_complement(background, words)
            rows = []
            keys = []
            for side, eq in enumerate(self.cuts):
                for u, v in product(range(10), repeat=2):
                    if weight(PAIRS[EXTRA[u]] + PAIRS[EXTRA[v]]) != char:
                        continue
                    keys.append((side, u, v))
                    entries = []
                    for w in words:
                        i, j = cut(w, side)
                        entries.append(eq.left[u, i] * eq.right[j, v])
                    rows.append(entries)
            linear = sp.Matrix(len(keys), len(words), sum(rows, [])) * normal
            pivot = list(linear.rref()[1])
            free = [i for i in range(normal.cols) if i not in pivot]
            independent = list(linear[:, pivot].T.rref()[1])
            self.blocks.append((words, normal, keys, linear, free, pivot, independent))

    def run(self, choose_free):
        tensor = dict(self.background)
        forced = 0
        for words, normal, keys, linear, free, pivot, independent in self.blocks:
            residual = [eq.residual(tensor) for eq in self.cuts]
            rhs = -sp.Matrix(len(keys), 1, [residual[s][u, v] for s, u, v in keys])
            nonlinear = any(rhs)
            vals = sp.zeros(normal.cols, 1)
            for i in free:
                vals[i] = choose_free(words, normal, i)
            rhs -= linear*vals
            chosen = linear[:, pivot]
            targets = rhs.extract(independent, [0])
            solution = chosen[independent, :].inv() * targets
            if chosen*solution != rhs:
                return None, weight(words[0]), forced
            for j, i in enumerate(pivot):
                vals[i] = solution[j]
            forced += int(nonlinear)
            tensor.update({w: v for w, v in zip(words, normal*vals) if v})
        assert all(eq.residual(tensor) == sp.zeros(10) for eq in self.cuts)
        return tensor, None, forced


def verify_recursion():
    # A rank-six central extension of the triple-face quantum background.
    symbolic, corrections = basis_tensors()
    base = {w: sp.sympify(v).subs({'Q': 2}) for w, v in symbolic.items()}
    target = {}
    cocharacter = (6, -2, 1, 1)
    for i, tensor in enumerate(corrections):
        if dot(weight(next(iter(tensor))), cocharacter) > 0:
            target.update({w: (i+1)*v for w, v in tensor.items()})
    recursion = WeightRecursion(base, cocharacter)
    def recover(words, normal, i):
        values = sp.Matrix([target.get(w, 0) for w in words])
        return normal.gauss_jordan_solve(values)[0][i]
    result, obstruction, _ = recursion.run(recover)
    assert obstruction is None and result == {**base, **target}

    # Arbitrary tangent coordinates at the commutative point can be obstructed.
    base = {w: (-1)**sum(w[i] > w[j] for i, j in combinations(AXES, 2))
            for w in permutations(AXES)}
    recursion = WeightRecursion(base, CHAMBERS[1])
    _, obstruction, _ = recursion.run(lambda words, normal, i: i+1)
    assert obstruction == (0, -1, 1, 0)
    selected = {((0, -1, 1, 0), 7), ((0, 0, -1, 1), 8)}
    result, obstruction2, forced = recursion.run(
        lambda words, normal, i: int((weight(words[0]), i) in selected))
    assert result and obstruction2 is None and forced == 1
    return obstruction


if __name__ == '__main__':
    verify_chambers()
    obstruction = verify_recursion()
    print('PASS: 336 chambers, 72 maximal supports, exact reconstruction, commutative obstruction', obstruction)
