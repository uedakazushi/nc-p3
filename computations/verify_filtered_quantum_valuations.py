"""First Koszul determinant in every ordering of a positive quantum deformation."""
from itertools import combinations, permutations, product

import sympy as sp

LETTERS = tuple(range(4))
CHAMBERS = ((0,1,4,6), (0,3,5,9), (0,6,8,9))
EXPONENTS = sp.symbols('c:6')


def verify(weights, order):
    rank = {j:i for i,j in enumerate(order)}
    oriented = dict(zip(combinations(LETTERS,2), EXPONENTS))
    def value(j,k):
        return oriented[j,k] if j<k else -oriented[k,j]
    def word_order(word):
        return sum(value(j,k) for h,j in enumerate(word)
                   for k in word[h+1:] if rank[j]<rank[k])
    # Only words of weight strictly above the volume allow arbitrary coefficients.
    # Equal weights have precisely the background support.
    triples = [tuple(j for j in reversed(order) if j!=k) for k in LETTERS]
    allowed = [[weights[a]>weights[b] or a==b for b in LETTERS] for a in LETTERS]
    surviving = [p for p in permutations(LETTERS) if all(allowed[i][p[i]] for i in LETTERS)]
    assert surviving == [LETTERS]
    degree = sum(word_order((a,)+triples[a]) for a in LETTERS)
    predicted = sum(value(j,k) for i,j in enumerate(order) for k in order[i+1:])
    assert sp.expand(degree-predicted) == 0
    assert sp.expand(word_order(tuple(reversed(order)))) == 0
    # Reversing the order and exchanging one adjacent pair recover all pairwise valuations.
    for h in range(3):
        other = list(order)
        other[h],other[h+1] = other[h+1],other[h]
        new = sum(value(j,k) for i,j in enumerate(other) for k in other[i+1:])
        assert sp.expand(predicted-new-2*value(order[h],order[h+1])) == 0


if __name__ == '__main__':
    for weights,order in product(CHAMBERS+((0,0,1,1), (0,0,0,0)), permutations(LETTERS)):
        verify(weights,order)
    print('All 24 monic determinants per chamber, including tied weights and adjacent transpositions.')
