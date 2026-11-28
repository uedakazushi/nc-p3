"""Remove the twenty rank-one contributions from a binary resultant on the kernel target."""

from itertools import combinations_with_replacement

from flint import nmod_mat, nmod_poly

from verify_configuration_intersections import (
    PAIRS, change_of_square_coordinates, exchange, modular_points, partner, rational_example,
)


def macaulay(quadrics, prime):
    words = list(combinations_with_replacement(range(4), 5))
    indices = {word: i for i, word in enumerate(words)}
    rows, mixed = [], []
    for index, word in enumerate(words):
        squares = [i for i in range(4) if word.count(i) >= 2]
        if len(squares) > 1:
            mixed.append(index)
        i = squares[0]
        rest = list(word)
        rest.remove(i)
        rest.remove(i)
        row = [0]*len(words)
        for pair, coefficient in zip(PAIRS, quadrics[i]):
            row[indices[tuple(sorted([*rest, *pair]))]] = coefficient
        rows.append(row)
    denominator = int(nmod_mat([[rows[i][j] for j in mixed] for i in mixed], prime).det())
    if denominator:
        return int(nmod_mat(rows, prime).det())*pow(denominator, -1, prime) % prime
    return None


def residual(first, third, prime):
    change = change_of_square_coordinates(first, third, prime)
    start, direction = [2, 17, 11, 5], [31, 7, 19, 29]
    def at(target):
        return [[sum(target[j]*int(change[k, PAIRS.index(tuple(sorted((i,j))))]) for j in range(4))
                 % prime for k in range(10)] for i in range(4)]
    a, b = at(start), at(direction)
    samples = []
    for t in range(prime):
        value = macaulay([[(s+t*d) % prime for s, d in zip(v,w)] for v,w in zip(a,b)], prime)
        if value is not None:
            samples.append((t,value))
        if len(samples) == 35:
            break
    assert len(samples) == 35
    matrix = nmod_mat([[pow(t, i, prime) for i in range(33)] for t,_ in samples[:33]], prime)
    coeff = matrix.solve(nmod_mat([[v] for _,v in samples[:33]], prime))
    polynomial = nmod_poly([int(coeff[i,0]) for i in range(33)], prime)
    assert polynomial.degree() == 32 and all(polynomial(t) == v for t,v in samples[33:])
    vertical = nmod_poly([1], prime)
    for point in third:
        factor = [sum(int(point[i])*v[i] for i in range(4)) % prime for v in (start,direction)]
        vertical *= nmod_poly(factor, prime)
    assert vertical.degree() == 20
    quotient, remainder = divmod(polynomial, vertical)
    assert not remainder and quotient.degree() == 12
    return quotient


if __name__ == '__main__':
    prime = 2657
    tensor, _, _ = rational_example()
    a = modular_points(tensor, prime)
    b = partner(a, tensor, False, prime)
    c = partner(b, exchange(tensor), True, prime)
    d = partner(c, tensor.T, False, prime)
    for first, third in ((a,c), (b,d), (c,a), (d,b)):
        q = residual(first, third, prime)
        assert q.gcd(q.derivative()).degree() == 0
    print("All four degree-32 resultants have square-free residual degree twelve after the twenty linear factors.")
