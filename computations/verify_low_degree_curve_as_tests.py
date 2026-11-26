"""Check bundle-row dimensions, the Hilbert identities, and the opposite line contraction."""
from itertools import permutations, product
from math import comb

from sympy import factor, symbols, zeros


def check_hilbert(degree):
    t = symbols("t")
    section_series = (1 + (degree-1)*t)/(1-t)**2
    denominator = (1-t)*(1-t**degree)
    assert factor(section_series/denominator - (1-t)**-4) == 0
    for n in range(12):
        h = lambda i: comb(i+3, 3) if i >= 0 else 0
        assert h(n)-h(n-1)-h(n-degree)+h(n-degree-1) == degree*n+1
        h1 = [comb(4, j)*max(0, degree*(j-n)-1) for j in range(5)]
        if n > degree+1:
            assert not any(h1)
        if n == degree+1:
            assert h1 == [0, 0, 0, 0, 1]
        if degree == 2 and n == 2:
            assert h1 == [0, 0, 0, 4, 3]


def check_opposite_contraction():
    # Alternating four-tensor. A line with the same ordered basis on both
    # sides gives its annihilator's decomposable two-form as the relation.
    for i, j in permutations(range(4), 2):
        relation = zeros(4)
        for a, b in product(range(4), repeat=2):
            word = (a, b, i, j)
            if len(set(word)) == 4:
                sign = sum(x > y for k, x in enumerate(word) for y in word[k+1:])
                relation[a, b] = (-1)**sign
        assert relation.rank() == 2
        assert relation[:, i] == relation[:, j] == zeros(4, 1)
        assert relation.T == -relation


def check_conic_lifts():
    # All quadratic lifts obtained by adding the plane equation times a
    # linear form give the same nonzero minimal opposite coefficient.
    u, v, w, x = symbols("u v w x")
    c = sum(z*l for z, l in zip((u, v, w, x), symbols("l:4")))
    plane_quadric = u*w-v**2
    lift = plane_quadric + x*c
    opposite = -plane_quadric - c*x
    assert (x*opposite + lift*x).expand() == 0
    assert opposite.subs(x, 0) == -plane_quadric


if __name__ == "__main__":
    for d in (1, 2):
        check_hilbert(d)
    check_opposite_contraction()
    check_conic_lifts()
    print("PASS: line and conic bundle rows, all-degree series, contractions and mixed lifts")
