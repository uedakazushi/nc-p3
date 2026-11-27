"""Exact cubic overlaps and line evaluation for a filtered four-letter example."""
from itertools import combinations, product

import sympy as sp

a, b, c, x, y = sp.symbols("a b c x y")
rules = {(j, i): {(i, j): 1} for i, j in combinations(range(4), 2)}
rules[2, 1][0, 0] = a
rules[3, 1][0, 0] = b
rules[3, 2][0, 0] = c


def order(word):
    return sum(letter >= 2 for letter in word), word


def reduce_word(word, from_right):
    descents = [i for i in range(len(word) - 1) if word[i:i+2] in rules]
    if not descents:
        return {word: 1}
    i = descents[-1 if from_right else 0]
    answer = {}
    for pair, coeff in rules[word[i:i+2]].items():
        child = word[:i] + pair + word[i+2:]
        assert order(child) < order(word)
        for normal, value in reduce_word(child, from_right).items():
            answer[normal] = sp.expand(answer.get(normal, 0) + coeff * value)
    return {key: value for key, value in answer.items() if value != 0}


words = list(product(range(4), repeat=3))
normal = set()
for word in words:
    left = reduce_word(word, False)
    assert left == reduce_word(word, True)
    normal.update(left)
assert len(normal) == 20
assert sum(word[0] > word[1] > word[2] for word in words) == 4

relations = []
for pair, tail in rules.items():
    coeff = sp.zeros(4)
    coeff[pair] = 1
    for lower, value in tail.items():
        coeff[lower] -= value
    relations.append(coeff)
for matrices in (relations, [matrix.T for matrix in relations]):
    evaluation = sp.Matrix.hstack(*(matrix.T * sp.Matrix([0, 0, x, y]) for matrix in matrices))
    assert evaluation.rank() == 3
    minors = [evaluation.extract(rows, cols).det()
              for rows in combinations(range(4), 3) for cols in combinations(range(6), 3)]
    assert sp.gcd_list(minors) == 1
    mixed = sp.Matrix([list(matrix[2:, :2]) for matrix in matrices[:5]])
    assert mixed.rank() == 4
assert relations[0][:2, :2].det() == 1
print("PASS: both line evaluations, both mixed projections and all cubic overlaps over Q[a,b,c]")
