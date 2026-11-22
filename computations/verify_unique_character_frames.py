"""Recover the first-edge frame and check the surviving monomial changes exactly."""

from itertools import product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_decomp

from verify_klein_elliptic_completion import gram_twist


def first_branch(w):
    thick = [
        v.reshape(4, 4)
        for v in w.columnspace()
        if not any(v[4 * i + j] for i in range(2) for j in range(2, 4))
    ]
    low = sp.Matrix.hstack(*[m[:2, :2].reshape(4, 1) for m in thick])
    high = sp.Matrix.hstack(*[m[2:, 2:].reshape(4, 1) for m in thick])
    t, s = sp.symbols("t s")
    matrix = (high * low.inv() * sp.Matrix([1, s, t, t * s])).reshape(2, 2)
    form = sp.Poly(sp.discriminant(matrix.det(), s), t)
    assert form.degree() == 4 and form.discriminant() != 0
    assert form.nth(1) == form.nth(3) == 0 and form.nth(0) ** 2 != form.nth(4) ** 2
    return sp.factor(form.nth(0) / form.nth(4))


def changes():
    for e, a, b, c in product(range(2), range(4), range(2), range(2)):
        shift = (a, a ^ b, a ^ c, a ^ b ^ c)
        yield tuple(tuple(i ^ u ^ (e * (i >> 1)) for i in range(4)) for u in shift)


def torus_relations(support):
    weights = sp.Matrix(
        [
            [int(word[p] == i) for p, i in product(range(4), repeat=2)]
            for word in support
        ]
    )
    smith, moves, _ = smith_normal_decomp(weights, domain=sp.ZZ)
    rank = sum(smith[j, j] != 0 for j in range(min(smith.shape)))
    assert [smith[j, j] for j in range(rank)] == [1] * 11 + [2, 2]
    assert moves[rank:, :] * weights == sp.zeros(len(support) - rank, 16)
    return [
        [(j, int(v)) for j, v in enumerate(row) if v]
        for row in moves[rank:, :].tolist()
    ]


def stabilizers(w):
    entries = {
        word: w[4 * word[0] + word[1], 4 * word[2] + word[3]]
        for word in product(range(4), repeat=4)
    }
    support = [word for word, v in entries.items() if v]
    relations = torus_relations(support)
    survivors = []
    for perm in changes():
        moved = [entries[tuple(q[i] for q, i in zip(perm, word))] for word in support]
        if all(moved):
            ratios = [a / entries[word] for a, word in zip(moved, support)]
            if all(sp.prod(ratios[j] ** e for j, e in row) == 1 for row in relations):
                survivors.append(perm)
    return survivors


def check():
    r = sp.Rational
    samples = [
        (
            sp.Matrix([[r(36, 17), 0], [9, 4]]),
            sp.Matrix([[r(36, 19), r(36, 37)], [r(36, 5), -12]]),
            3,
            2,
            r(3195075625, 15968016),
        ),
        (
            sp.Matrix([[1, 1], [1, 2]]),
            sp.Matrix([[-r(13, 14), r(1, 13)], [-14, -r(16, 11)]]),
            -2,
            -2,
            r(9604, 184041),
        ),
    ]
    assert len(set(changes())) == 32
    for a, b, z, l, expected in samples:
        w, _, _ = gram_twist(a, b, sp.Integer(z), sp.Integer(l))
        assert first_branch(w) == expected
        assert stabilizers(w) == [(tuple(range(4)),) * 4]
        print(
            f"z={z}: unique compatible frame, branch ratio {expected}, one permutation in the torus stabilizer."
        )


if __name__ == "__main__":
    check()
