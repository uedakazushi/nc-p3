"""The noncommutative Ore overlap at a pair of double quantum faces."""

from collections import defaultdict
from itertools import combinations, permutations

import sympy as sp

a, A, B, C, D, H, u, v, s, t = sp.symbols("a A B C D H u v s t")
a = sp.symbols("a", nonzero=True)
x, y, z = (0,), (1,), (2,)


def add(values, word, coeff):
    values[word] += coeff


def normal(word, coeff, result):
    for i, (first, second) in enumerate(zip(word, word[1:])):
        if first <= second:
            continue
        left, right = word[:i], word[i+2:]
        if (first, second) == (1, 0):
            normal(left + x+y + right, coeff, result)
            normal(left + z+z + right, coeff*D, result)
        else:
            scalar = a if second == 0 else 1/a
            normal(left+(second, first)+right, coeff*scalar, result)
        return
    add(result, word, coeff)


def verify_ore(gamma):
    delta = ({x+z:u, y+z:v}, {x+z:s, y+z:t}, {x+x:A, x+y:B, y+y:C, z+z:H})
    sigma = (1/a, a, gamma)
    relations = ({z+x:1, x+z:-a}, {z+y:1, y+z:-1/a}, {y+x:1, x+y:-1, z+z:-D})
    checks = []
    for relation in relations:
        coeffs = defaultdict(int)
        for (i, j), k in relation.items():
            for word, coeff in delta[i].items():
                normal(word+(j,), k*coeff, coeffs)
            for word, coeff in delta[j].items():
                normal((i,)+word, k*sigma[i]*coeff, coeffs)
        checks.append({word:sp.factor(coeff) for word, coeff in coeffs.items() if sp.factor(coeff)})
    subs = {s:D*A*(1+gamma*a**2)/(a-1/a), v:D*C*(1+gamma/a**2)/(a-1/a)}
    if gamma == 1:
        subs.update({H:-D*B/(a**2-1), u:0, t:2*D*B/(a-1/a)})
    else:
        subs.update({H:0, u:D*B/(2*a), t:-D*B/(2*a)})
    residual = {(i, word):sp.factor(coeff.subs(subs))
                for i, check in enumerate(checks) for word, coeff in check.items()}
    residual = {term:coeff for term, coeff in residual.items() if coeff}
    expected = {(2, z*3):2*D**2*B*(a*a+1)/(a*a-1)} if gamma == 1 else {}
    assert set(residual) == set(expected)
    assert all(sp.factor(residual[k]-v) == 0 for k, v in expected.items())


def verify_supports():
    # Residue exponents mod 24, with product signs represented by 0 and 12.
    weights = [tuple(i if n == 0 else 2-i if n == 1 else -1 for n in range(4)) for i in range(3)]
    weights += [w[2:]+w[:2] for w in weights]
    for epsilon in (0, 12):
        faces = (5, 5, (epsilon-5) % 24, (epsilon-5) % 24)
        group = []
        for perm in permutations(range(4)):
            sign = (-1)**sum(perm[i]>perm[j] for i, j in combinations(range(4), 2))
            if all((faces[perm[i]]-sign*faces[i]) % 24 == 0 for i in range(4)):
                group.append(perm)
        assert len(group) == (4 if epsilon == 0 else 2)
        supports = []
        for side in (0, 3):
            for j in (0, 2):
                mixed = {side, side+1, side+2, (3-side)+j}
                if epsilon == 0:
                    supports.append(frozenset(mixed - {side+1}))
                    supports.append(frozenset({side, side+1, side+2}))
                else:
                    supports.append(frozenset(mixed))
        maximal = set(supports)
        def orbit(support):
            return {frozenset(weights.index(tuple(weights[i][j] for j in p)) for i in support) for p in group}
        assert set.union(*(orbit(s) for s in maximal)) == maximal
        reps = {min(tuple(sorted(s)) for s in orbit(supp)) for supp in maximal}
        dimensions = sorted(len(r) - sp.Matrix([weights[j] for j in r]).rank() for r in reps)
        assert dimensions == ([0, 1] if epsilon == 0 else [1, 1])


if __name__ == '__main__':
    verify_ore(1)
    verify_ore(-1)
    verify_supports()
    print('PASS: all ternary Ore overlaps, the cubic obstruction, finite permutations and maximal support dimensions')
