"""Certify the order-four graded root and its central cuts over Q(i) and F_181."""
from itertools import product

from flint import fmpq_mat, nmod_mat
from sympy.polys.domains import QQ_I

from verify_elliptic_quintic_overlaps import IMAGINARY, PRIME, torsion_relations

ZERO, ONE, I = QQ_I.zero, QQ_I.one, QQ_I(0, 1)


def complex_rref(rows):
    """Use a native rational RREF with real and imaginary coordinates interleaved."""
    height, width = len(rows), len(rows[0])
    matrix = fmpq_mat(2*height, 2*width)
    for i, row in enumerate(rows):
        for j, z in enumerate(row):
            if z:
                matrix[2*i, 2*j] = matrix[2*i+1, 2*j+1] = z.x
                matrix[2*i, 2*j+1], matrix[2*i+1, 2*j] = -z.y, z.y
    reduced, rank = matrix.rref()
    assert rank % 2 == 0
    return [[QQ_I(reduced[2*i, 2*j], reduced[2*i+1, 2*j]) for j in range(width)]
            for i in range(rank//2)]


def projection(relations):
    """Columns give a basis for the kernel, normalized on the free coordinates."""
    reduced = complex_rref(relations)
    pivots = {next(i for i, c in enumerate(row) if c): row for row in reduced}
    free = [i for i in range(len(relations[0])) if i not in pivots]
    return [[-pivots[i][j] if i in pivots else QQ_I(int(i == j)) for j in free]
            for i in range(len(relations[0]))], free


def root_relations():
    relations = torsion_relations(1)
    diagonal = [ONE, -ONE, I, -I]
    return [[QQ_I.from_sympy(relations[k, h])*diagonal[k % 4] for k in range(16)] for h in range(6)]


def frobenius(relations):
    normals, _ = projection(relations)
    frame = relations
    for degree, dimension in [(3, 4), (4, 1)]:
        conditions = []
        for prefix, normal in product(range(4**(degree-2)), range(10)):
            conditions.append([sum((vector[4*prefix+a]*normals[4*a+b][normal] for a in range(4)), ZERO)
                               for vector in frame for b in range(4)])
        kernel, free = projection(conditions)
        assert len(free) == dimension
        frame = [[sum((vector[word//4]*kernel[4*j+word%4][h] for j, vector in enumerate(frame)), ZERO)
                  for word in range(4**degree)] for h in range(dimension)]
    top, = frame
    assert all(c == top[(i % 64)*4+i//64] for i, c in enumerate(top))
    for cut, rank in [(1, 4), (2, 6), (3, 4)]:
        width = 4**(4-cut)
        assert len(complex_rref([top[i:i+width] for i in range(0, 256, width)])) == rank


def exact_paths(relations):
    words, arrows = [[()]], []
    for degree in range(1, 6):
        if degree == 1:
            values = [[QQ_I(int(i == j)) for j in range(4)] for i in range(4)]
            free = range(4)
        else:
            rows = []
            for i, rule in product(range(len(words[-2])), relations):
                coefficients = [ZERO]*(4*len(words[-1]))
                for b, c in product(range(4), repeat=2):
                    scalar = rule[4*b+c]
                    if scalar:
                        for j, coeff in enumerate(arrows[-1][4*i+b]):
                            if coeff:
                                coefficients[4*j+c] += scalar*coeff
                rows.append(coefficients)
            values, free = projection(rows)
        words.append([words[-1][i//4]+(i%4,) for i in free])
        arrows.append(values)
    return words, arrows


def word_value(arrows, word):
    result = [ONE]
    for degree, letter in enumerate(word):
        matrix = arrows[degree]
        result = [sum((c*matrix[4*i+letter][j] for i, c in enumerate(result)
                       if c and matrix[4*i+letter][j]), ZERO) for j in range(len(matrix[0]))]
    return result


def quartics(relations):
    words, arrows = exact_paths(relations)
    assert [len(group) for group in words] == [1, 4, 10, 20, 35, 56]
    left = [[word_value(arrows, (b,)+w) for b in range(4)] for w in words[4]]
    result = []
    for scalar, dimension in [(ONE, 2), (-ONE, 1), (I, 1), (-I, 1)]:
        system = [[arrows[4][4*i+b][j]-scalar*left[i][b][j] for i in range(35)]
                  for b in range(4) for j in range(56)]
        kernel, free = projection(system)
        assert len(free) == dimension
        for h in range(dimension):
            result.append((scalar, {w: kernel[i][h] for i, w in enumerate(words[4]) if kernel[i][h]}))
    return result


def residue(z):
    def reduce_rational(value):
        return int(value.numerator)*pow(int(value.denominator), -1, PRIME) % PRIME
    return (reduce_rational(z.x)+IMAGINARY*reduce_rational(z.y)) % PRIME


def folded_cuts(quartic_basis):
    support = sorted(set().union(*(coefficients for _, coefficients in quartic_basis)))
    return [[{w: sum(residue(scalar)**vertex * k**g * residue(coefficients.get(w, ZERO))
                     for k, (scalar, coefficients) in enumerate(quartic_basis, 1)) % PRIME
              for w in support} for g in range(4)] for vertex in range(4)]


def cut_dimensions(relations, cuts, start):
    dimensions, arrows = [1], []
    for degree in range(1, 14):
        width = 4*dimensions[-1]
        rows = []
        for length, rules in [(2, relations), (4, cuts[(start+degree-4) % 4])]:
            if length > degree:
                continue
            for rule in rules:
                components = [nmod_mat(dimensions[-length], dimensions[-1], PRIME) for _ in range(4)]
                for word, scalar in rule.items():
                    if scalar:
                        value = arrows[degree-length][word[0]]
                        for k, b in enumerate(word[1:-1], 1):
                            value = value*arrows[degree-length+k][b]
                        components[word[-1]] += value*scalar
                rows.extend([[components[b][i, j] for j in range(dimensions[-1]) for b in range(4)]
                             for i in range(dimensions[-length])])
        matrix = nmod_mat(rows, PRIME) if rows else nmod_mat(0, width, PRIME)
        kernel, rank = matrix.nullspace()
        arrows.append([nmod_mat(dimensions[-1], rank,
                                [kernel[4*i+b, j] for i in range(dimensions[-1]) for j in range(rank)], PRIME)
                       for b in range(4)])
        dimensions.append(rank)
    return dimensions


def main():
    relations = root_relations()
    frobenius(relations)
    cuts = folded_cuts(quartics(relations))
    rules = [{(j//4, j%4): residue(c) for j, c in enumerate(row) if c} for row in relations]
    expected = [1, 4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1, 0]
    for start in range(4):
        assert cut_dimensions(rules, cuts, start) == expected
    print("PASS: exact Frobenius pairings, five exact quartics, and four central reductions")


if __name__ == "__main__":
    main()
