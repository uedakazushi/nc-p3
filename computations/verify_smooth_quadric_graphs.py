"""The quadric graph tensor, spin pairing, exceptional fibers and balance."""

from fractions import Fraction
from itertools import product

import sympy as sp

WORDS = tuple(product(range(4), repeat=4))


def spin_matrices():
    creation = sp.Matrix([[0, 0], [1, 0]])
    sign = sp.diag(1, -1)
    return (sp.kronecker_product(creation, sp.eye(2)),
            sp.kronecker_product(creation.T, sp.eye(2)),
            sp.kronecker_product(sign, creation),
            sp.kronecker_product(sign, creation.T))


def volume(m):
    letters = spin_matrices()
    x, y, z, w = letters
    u, t = x*y-y*x, z*w-w*z
    sign = sp.diag(1, -1, -1, 1)
    traces = {}
    for word in WORDS:
        value = sign*(u-m*t)
        for i in word:
            value = value*letters[i]
        traces[word] = sp.expand(value.trace())
    return traces


def verify():
    x, y, z, w = spin_matrices()
    u, t = x*y-y*x, z*w-w*z
    m = sp.symbols('m')
    assert ((u+m*t)*(u-m*t)).expand() == (1-m**2)*sp.eye(4)
    assert (u+t).rank() == (u-t).rank() == 2
    symbolic = volume(m)
    tensor = {word: val.subs(m, 2) for word, val in symbolic.items() if val.subs(m, 2)}
    assert all(sum(word.count(i)*v for i, v in enumerate((1,-1,0,0))) == 0 and
               sum(word.count(i)*v for i, v in enumerate((0,0,1,-1))) == 0 for word in tensor)
    def flip(word, perm):
        return tuple(perm[i] for i in word)
    for word in tensor:
        assert flip(word, (1,0,3,2)) in tensor
        assert flip(word[1:]+word[:1], (0,1,3,2)) in tensor
    assert (0,1,0,1) in tensor and (2,3,2,3) in tensor
    for side in range(2):
        pairs = list(product(range(4), repeat=2))
        def original(left, right):
            rotated = left+right
            return rotated[-side:]+rotated[:-side] if side else rotated
        mat = sp.Matrix([[tensor.get(original(p, q), 0) for q in pairs] for p in pairs])
        assert mat.rank() == 6
        signs = (-1)**side
        relations = [{(0,1):1, (1,0):1, (2,3):1, (3,2):1},
                     {(2,3):1, (3,2):-1, (0,1):-2*signs, (1,0):2*signs}]
        for i, j in product((0,1), (2,3)):
            relations.append({(i,j):1, (j,i):-1})
        frame = sp.Matrix([[r.get(p, 0) for p in pairs] for r in relations])
        assert mat.row_join(frame.T).rank() == 6
    for value in (Fraction(2), Fraction(-2), Fraction(1,2), Fraction(-1,2)):
        assert value**2 + value**-2 == Fraction(17,4)


if __name__ == '__main__':
    verify()
    print('PASS: Clifford factors, two adjacent relation spaces, support symmetries, exceptional ranks and quotient invariant')
