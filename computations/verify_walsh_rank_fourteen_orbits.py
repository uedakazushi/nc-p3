"""Intrinsic four-line frames, unrestricted isomorphisms and periods at Walsh rank fourteen."""

from functools import cache
from itertools import product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_decomp

from verify_rank_fifteen_parameter_orbits import pattern, template_group


LABELS = range(4)
R, Q, N = sp.symbols("r q n")
ATOMS = (1, R, Q, N)


def coefficients(negative):
    p = -R if negative else R
    m = l = sp.I * R
    return pattern(1, 0, (p+m)/2, (p-m)/2, (N-Q)/2, (N+Q)/2, (l-R)/2, (l+R)/2)


def verify_marked_orbits():
    group = [(p, e) for p, e in template_group() if p[1] == 1]
    for negative in (False, True):
        gamma = coefficients(negative)
        rows = set()
        for perm, exponents in group:
            row = [gamma[j] * sp.I**e for j, e in zip(perm, exponents, strict=True)]
            assert row[8] in (-1, 1)
            rows.add(tuple(sp.expand(v / row[8]) for v in row))
        assert len(rows) == 64


def framed_tensor(negative):
    gamma = coefficients(negative)
    factors = {sp.expand(sp.I**k * v): (k, j) for k in LABELS for j, v in enumerate(ATOMS)}
    entries = {}
    for word in product(LABELS, repeat=4):
        i, j, k, t = word
        if i ^ k != j ^ t:
            continue
        d, g = (i ^ k) & 1, (i ^ k) >> 1
        value = sp.expand(sp.S.Half * sum(
            gamma[h + 2*g + 4*(d + 2*b)]
            * (-1) ** (((i ^ j) & 1)*h + b*(k >> 1) + (b ^ d)*(t >> 1))
            for h, b in product(range(2), repeat=2)
        ))
        if value != 0:
            entries[word] = factors[value]
    return entries


def frame_permutations():
    for flip, first, second in product(range(2), LABELS, LABELS):
        def perm(bits):
            return tuple(i ^ flip ^ (((bits >> (i & 1)) & 1) * 2) for i in LABELS)
        p, q = perm(first), perm(second)
        yield p, q, p, q


@cache
def relations(support):
    weights = sp.Matrix([[int(w[p] == i) for p, i in product(LABELS, repeat=2)] for w in support])
    smith, moves, _ = smith_normal_decomp(weights, domain=sp.ZZ)
    rank = sum(smith[i, i] != 0 for i in range(min(smith.shape)))
    kernel = moves[rank:, :]
    assert kernel * weights == sp.zeros(len(support) - rank, 16)
    return [[(j, int(v)) for j, v in enumerate(row) if v] for row in kernel.tolist()]


def verify_comparison(negative, shift, zero_atom=None):
    source = framed_tensor(False)
    target = framed_tensor(negative)
    if zero_atom is not None:
        source = {w: z for w, z in source.items() if z[1] != zero_atom}
        target = {w: z for w, z in target.items() if z[1] != zero_atom}
    support = tuple(sorted(source))
    counts = [0, 0, 0]
    for perm in frame_permutations():
        ratios = []
        for word in support:
            moved = tuple(p[i] for p, i in zip(perm, word, strict=True))
            moved = moved[shift:] + moved[:shift]
            if moved not in target:
                break
            ratios.append((source[word], target[moved]))
        else:
            counts[0] += 1
            checks = set()
            for relation in relations(support):
                powers = [0] * 4
                for j, exponent in relation:
                    (p, a), (q, b) = ratios[j]
                    powers[0] += exponent * (q - p)
                    if a:
                        powers[a] -= exponent
                    if b:
                        powers[b] += exponent
                powers[0] %= 4
                checks.add(tuple(powers))
            if all(not any(row) for row in checks):
                counts[1] += 1
            elif any(k and not (r or q or n) for k, r, q, n in checks):
                counts[2] += 1
    return tuple(counts)


if __name__ == "__main__":
    verify_marked_orbits()
    assert verify_comparison(True, 0) == (8, 0, 8)
    assert verify_comparison(False, 1) == (8, 0, 8)
    assert verify_comparison(True, 1) == (8, 4, 4)
    assert verify_comparison(False, 2) == (8, 4, 4)
    for missing in (2, 3):
        assert verify_comparison(True, 1, missing)[1] > 0
    print("PASS: rank-fourteen sheets, intrinsic monomial equivalence and shift")
