"""Check degree-16 Koszul exactness at 181 using the five character spaces."""
from functools import cache
from itertools import product
from math import comb

from flint import nmod_mat

from verify_elliptic_quintic_overlaps import COSETS, PRIME, add, coordinates, nullvectors


def matrix(rows, width):
    return nmod_mat(rows, PRIME) if rows else nmod_mat(0, width, PRIME)


def concatenate(blocks):
    width = sum(m.ncols() for m in blocks)
    height = blocks[0].nrows()
    return nmod_mat(height, width, [m[i, j] for i in range(height)
                                    for m in blocks for j in range(m.ncols())], PRIME)


def stack(blocks, width):
    return nmod_mat(sum(m.nrows() for m in blocks), width,
                    [x for m in blocks for x in m.entries()], PRIME)


@cache
def frame(position, length):
    """All tensors satisfying each adjacent evaluation, grouped by weight."""
    words = list(product(range(1, 5), repeat=length))
    vectors = []
    for weight in range(5):
        block = [w for w in words if sum(w) % 5 == weight]
        rows = []
        for k in range(length-1):
            shift = [(1, 0), (PRIME-1, 0)][(position+k) % 2]
            outside = sorted({w[:k]+w[k+2:] for w in block})
            for q in COSETS[:3]:
                u, v = coordinates(q), coordinates(add(q, shift))
                for fixed in outside:
                    rows.append([u[w[k]-1]*v[w[k+1]-1] % PRIME if w[:k]+w[k+2:] == fixed else 0
                                 for w in block])
        for coeffs in nullvectors(matrix(rows, len(block))):
            vectors.append((weight, {w: c for w, c in zip(block, coeffs) if c}))
    return vectors


def cubic_slices(position):
    quadrics = frame(1-position, 2)
    result = []
    for weight, cubic in frame(position, 3):
        coefficients = {}
        for b in range(1, 5):
            vectors = [v for h, v in quadrics if h == (weight-b) % 5]
            words = sorted(set().union(*vectors))
            basis = matrix([[v.get(w, 0) for w in words] for v in vectors], len(words))
            reduced, rank = basis.rref()
            pivots = [next(j for j in range(len(words)) if reduced[i, j]) for i in range(rank)]
            submatrix = matrix([[v.get(words[j], 0) for j in pivots] for v in vectors], rank)
            target = matrix([[cubic.get((b,)+words[j], 0) for j in pivots]], rank)
            values = target*submatrix.inv()
            assert (values*basis).entries() == [cubic.get((b,)+w, 0) for w in words]
            coefficients[b] = [int(x) for x in values.entries()]
        result.append((weight, coefficients))
    return result


def path_step(degree, dimensions, arrows):
    next_dimensions, step = [], {}
    for weight in range(5):
        widths = [dimensions[-1][(weight-b) % 5] for b in range(1, 5)]
        blocks = []
        if degree >= 2:
            for h, rule in frame((degree-2) % 2, 2):
                start = (weight-h) % 5
                parts = []
                for c in range(1, 5):
                    part = nmod_mat(dimensions[-2][start], widths[c-1], PRIME)
                    for (b, last), value in rule.items():
                        if last == c:
                            part += arrows[-1][start, b]*value
                    parts.append(part)
                blocks.append(concatenate(parts))
        kernel, dimension = stack(blocks, sum(widths)).nullspace()
        next_dimensions.append(dimension)
        offset = 0
        for b, width in enumerate(widths, 1):
            step[(weight-b) % 5, b] = nmod_mat(width, dimension,
                [kernel[i, j] for i in range(offset, offset+width) for j in range(dimension)], PRIME)
            offset += width
    return next_dimensions, step


def third_rank(degree, dimensions, arrows, slices):
    quadrics = frame((degree-2) % 2, 2)
    total = 0
    for weight in range(5):
        blocks = []
        for h, coefficients in slices:
            start = (weight-h) % 5
            parts, count = [], {}
            for k, _ in quadrics:
                j = count.get(k, 0)
                count[k] = j+1
                b = (h-k) % 5
                if b and coefficients[b][j]:
                    part = arrows[degree-3][start, b]*coefficients[b][j]
                else:
                    part = nmod_mat(dimensions[degree-3][start], dimensions[degree-2][(weight-k) % 5], PRIME)
                parts.append(part)
            blocks.append(concatenate(parts))
        width = sum(dimensions[degree-2][(weight-k) % 5] for k, _ in quadrics)
        total += stack(blocks, width).rank()
    return total


def main():
    dimensions, arrows = [[1, 0, 0, 0, 0]], []
    slices = [cubic_slices(i) for i in range(2)]
    for position in range(2):
        assert [len(frame(position, n)) for n in (2, 3, 4)] == [6, 4, 1]
        _, top = frame(position, 4)[0]
        for cut, expected_rank in [(1, 4), (2, 6), (3, 4)]:
            prefixes = list(product(range(1, 5), repeat=cut))
            suffixes = list(product(range(1, 5), repeat=4-cut))
            values = [[top.get(a+b, 0) for b in suffixes] for a in prefixes]
            assert matrix(values, len(suffixes)).rank() == expected_rank
    for degree in range(1, 17):
        sizes, step = path_step(degree, dimensions, arrows)
        dimensions.append(sizes)
        arrows.append(step)
        assert sum(sizes) == comb(degree+3, 3)
        previous = sum(dimensions[degree-4]) if degree >= 4 else 0
        if degree >= 4:
            for weight in range(5):
                action = concatenate([arrows[degree-4][weight, b] for b in range(1, 5)])
                assert action.rank() == dimensions[degree-4][weight]
        if degree >= 3:
            assert third_rank(degree, dimensions, arrows, slices[(degree-3) % 2]) == (
                4*sum(dimensions[degree-3])-previous)
    print("PASS: all Hilbert dimensions and Koszul ranks through degree 16")


if __name__ == "__main__":
    main()
