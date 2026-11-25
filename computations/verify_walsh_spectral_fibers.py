"""Recover the squared Walsh parameters from labeled Schmidt spectra, including ties."""

from itertools import combinations, permutations, product

import sympy as sp


PAIRS = ((0, 1), (1, 2), (2, 0))


def norm(value):
    return sp.simplify(sp.expand_complex(value * sp.conjugate(value)))


def phase(index, character):
    return (-1) ** (index & character).bit_count()


def labeled_spectra(matrix):
    rows = [sum(norm(z) for z in row) for row in matrix]
    cols = [sum(norm(matrix[i][j]) for i in range(3)) for j in range(3)]
    singles = [norm(sum(phase(i, h)*matrix[i][i] for i in range(3))) for h in range(4)]
    pairs = [(norm(matrix[i][j] + matrix[j][i]), norm(matrix[i][j] - matrix[j][i]))
             for i, j in PAIRS]
    return rows, cols, singles, pairs


def gram_coordinates(gram):
    """Real planar Gram realizations, up to a common complex phase."""
    anchor = next((i for i in range(3) if gram[i, i] > 0), None)
    if anchor is None:
        return [(0, 0, 0)] if gram == sp.zeros(3) else []
    x = [sp.simplify(gram[i, anchor] / sp.sqrt(gram[anchor, anchor])) for i in range(3)]
    residual = sp.simplify(gram - sp.Matrix(x)*sp.Matrix(x).T)
    if any(residual[i, i] < 0 for i in range(3)):
        return []
    other = next((i for i in range(3) if residual[i, i] > 0), None)
    if other is None:
        return [tuple(x)] if residual == sp.zeros(3) else []
    y = [sp.simplify(residual[i, other] / sp.sqrt(residual[other, other])) for i in range(3)]
    if sp.simplify(residual - sp.Matrix(y)*sp.Matrix(y).T) != sp.zeros(3):
        return []
    return [tuple(sp.simplify(a + sign*sp.I*b) for a, b in zip(x, y)) for sign in (-1, 1)]


def square_matrices(rows, cols, singles, pairs):
    """One labeling; inputs are exact real algebraic squared singular values."""
    product_norm = sp.simplify(pairs[0][0]*pairs[0][1])
    if any(sp.simplify(a*b - product_norm) != 0 for a, b in pairs):
        return
    sums = [sp.simplify((a+b)/2) for a, b in pairs]
    gram = sp.zeros(3)
    for i in range(3):
        gram[i, i] = sp.simplify((rows[i]+cols[i] - sum(s for ij, s in zip(PAIRS, sums) if i in ij))/2)
    if any(gram[i, i] < 0 for i in range(3)) or sp.simplify(sum(singles)-4*sp.trace(gram)) != 0:
        return
    for i, j in PAIRS:
        gram[i, j] = gram[j, i] = sp.simplify(sum(phase(i ^ j, h)*z for h, z in enumerate(singles))/8)
    for diag in gram_coordinates(gram):
        squares = [z*z for z in diag]
        differences = [sp.simplify(squares[3-i-j] - squares[i] - squares[j]) for i,j in PAIRS]
        first = next((i for i,d in enumerate(differences) if d != 0), None)
        if first is None:
            js = [sp.sqrt(product_norm)]
        else:
            d, s = differences[first], sums[first]
            a = sp.simplify(s*s*(norm(d) + product_norm - s*s))
            b = sp.simplify(norm(d)*product_norm - a)
            if a < 0 or b < 0:
                continue
            js = {sp.simplify((i*sp.sqrt(a) + sp.I*j*sp.sqrt(b))/sp.conjugate(d))
                  for i,j in product((-1, 1), repeat=2)}
        for signed_difference in js:
            matrix = sp.diag(*squares)
            for (i,j), d in zip(PAIRS, differences):
                matrix[i,j] = sp.simplify((d+signed_difference)/2)
                matrix[j,i] = sp.simplify((d-signed_difference)/2)
            yield matrix


def pair_slots(spectrum):
    """Label the four single slots and the three pairs with their doubled multiplicities."""
    # Eigenvalues include zeros. Sort only to give the multiset one representation.
    values = tuple(sorted(spectrum, key=sp.default_sort_key))
    for single in {tuple(values[i] for i in indices) for indices in combinations(range(16), 4)}:
        rest = list(values)
        for z in single:
            rest.remove(z)
        if any(rest[i] != rest[i+1] for i in range(0, 12, 2)):
            continue
        half = rest[::2]
        for pairs in set(permutations(half)):
            slots = tuple((pairs[i],pairs[i+1]) for i in (0,2,4))
            if any(sp.default_sort_key(a) > sp.default_sort_key(b) for a,b in slots):
                continue
            if len({sp.simplify(a*b) for a,b in slots}) != 1:
                continue
            for singles in set(permutations(single)):
                yield singles, slots



def adjacent_slots(spectrum):
    positive = sorted((v for v in spectrum if v != 0), key=sp.default_sort_key)
    if len(positive) > 6 or len(positive) % 2 or positive[::2] != positive[1::2]:
        return ()
    values = [v/2 for v in positive[::2]]
    return set(permutations(values + [sp.S.Zero]*(3-len(values))))


def normalization_candidates(first, second, opposite):
    """Candidates (signed square matrix, product of the six off-diagonal entries).

    Supply the three lists of sixteen exact squared singular values in a common
    scale. The resulting projective list can include equivalent labelings.
    """
    if any(len(values) != 16 for values in (first, second, opposite)):
        raise ValueError("each flattening needs its sixteen squared singular values")
    if not (sum(first) == sum(second) == sum(opposite)) or sum(first) == 0:
        return set()
    row_slots, col_slots = adjacent_slots(first), adjacent_slots(second)
    out = set()
    for singles, pairs in pair_slots(opposite):
        for rows, cols in product(row_slots, col_slots):
            for matrix in square_matrices(rows, cols, singles, pairs):
                # Squaring and choosing planar Gram coordinates can add spurious roots.
                if any(sp.simplify(abs(matrix[i,j]) + abs(matrix[j,i]) - (a+b)/2) != 0
                       for (i,j), (a,b) in zip(PAIRS, pairs)):
                    continue
                radical = sp.sqrt(sp.prod(matrix[i,j] for i,j in permutations(range(3),2)))
                for joint in (radical, -radical):
                    out.add((tuple(matrix), joint))
    return out

def verify():
    x, y, a, b = sp.symbols("x y a b", real=True)
    u, v = x+sp.I*y, a+sp.I*b
    assert sp.expand(norm(u+v)*norm(u-v) - ((norm(u)+norm(v))**2 - 4*(a*x+b*y)**2)) == 0
    matrices = [
        sp.Matrix([[0,1,sp.I],[sp.I,0,sp.I],[-1,1,0]]),
        sp.zeros(3),
    ]
    # Construct the second example from signed squares with a non-real circulation.
    d0,d1,d2,j = map(sp.sympify, (1,4,-1, 2+4*sp.I))
    for u,v in PAIRS:
        d = (d0,d1,d2)
        total = d[3-u-v] - d[u] - d[v]
        sign = phase(u ^ v, 3)
        matrices[1][u,v] = sp.sqrt(sign*(total+j)/2)
        matrices[1][v,u] = sp.sqrt(sign*(total-j)/2)
    matrices[1][0,0],matrices[1][1,1],matrices[1][2,2] = 1,2,sp.I
    for matrix in matrices:
        data = labeled_spectra(matrix.tolist())
        expected = sp.Matrix(3,3,lambda i,j: sp.simplify(phase(int(i^j),3)*matrix[i,j]**2))
        answers = list(square_matrices(*data))
        if not any(x == expected for x in answers):
            # With zero diagonal the common phase of the circulation is free.
            assert expected.diagonal() == sp.zeros(1,3)
            assert any(x*expected[0,1] == expected*x[0,1] for x in answers)
    adjacent = [sp.S.Zero]*10 + [sp.S(4)]*6
    repeated = [sp.S.Zero]*4 + [sp.S(2)]*12
    quantum = normalization_candidates(adjacent, adjacent, repeated)
    assert quantum and all(all(entries[k] == 0 for k in (0,4,8)) for entries,p in quantum)
    assert any(all(x == 0 for x in s) for s,_ in pair_slots(repeated))


if __name__ == "__main__":
    verify()
    print("PASS: inverse spectral identities, repeated spectra, zero diagonal and complex circulation")
