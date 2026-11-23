"""Uniform signed unit-residue covers and the complete F9 unit grid."""

from itertools import combinations, product

import sympy as sp
from flint import fq_default_ctx, fq_default_poly_ctx

from verify_rank_fifteen_unit_residues import verify_unit_point


def unit_cover():
    v, b, Z = sp.symbols("v b Z")
    D, d = v**3-v, b**3-b

    def vanishes(expression):
        return sp.Poly(expression, v, b, modulus=3).is_zero

    def radicand(a, sign):
        return (v-a-sign*b)*sp.prod(v-j for j in range(3) if j != a % 3)

    original = [(1-v*v)*(b+sign*2*v) for sign in (1, -1)]
    original += [2*v*(1+v)*(b+sign*(1-v)) for sign in (1, -1)]
    original += [2*v*(1-v)*(b+sign*(1+v)) for sign in (1, -1)]
    for a, sign in product(range(3), (1, -1)):
        f = radicand(a, sign)
        assert sum(any(vanishes(g-k*f) for k in (1, -1)) for g in original) == 1
        for epsilon, shift in product((1, -1), range(3)):
            assert vanishes(f.subs(v, epsilon*v+shift)
                            -epsilon*radicand(epsilon*(a-shift) % 3, epsilon*sign))
    for sign in (1, -1):
        roots = [radicand(a, sign) for a in range(3)]
        assert vanishes(sum(roots)-sign*b)
        assert vanishes(sum(first*second for first, second in combinations(roots, 2))-2*sign*b*D)
        assert vanishes(sp.prod(roots)-D**2*(D-sign*d))
        assert vanishes(roots[0]-roots[1]+sign*b*(2*v-1))

    r = sp.symbols("r:3")
    four_sums = [sum(x*y for x, y in zip(signs, r, strict=True))
                 for signs in product((1, -1), repeat=3) if sp.prod(signs) == 1]
    s1 = sum(x*x for x in r)
    s2 = sum(x*x*y*y for x, y in combinations(r, 2))
    assert sp.expand(sp.prod(Z-value for value in four_sums)
                     -(Z**4-2*s1*Z**2-8*sp.prod(r)*Z+s1*s1-4*s2)) == 0
    l, a, s = sp.symbols("l a s")
    A, B = s*a, s*(l-a)
    cubic = A*(A-s)**2+(B*B+(s**3-s)*A)*(B+s)
    normalized = s**2*(a*a-(l+1)*a)-l*l*(l+1)
    assert sp.Poly(cubic+s**3*normalized, l, a, s, modulus=3).is_zero
    Y = 2*a-l-1
    assert sp.Poly(s**2*Y**2-(l+1)*(s**2*(l+1)+l*l)-normalized,
                   l, a, s, modulus=3).is_zero
    X, t = sp.symbols("X t")
    assert sp.Poly((l+1)*(l+1+X*l*l)
                   -(X*t**3+(X+1)*t*t+X*t).subs(t, l+1), l, X, modulus=3).is_zero
    # c4=16(a^2-3), discriminant=16(a^2-4) for y^2=t^3+a*t^2+t.
    assert sp.Poly((X+1)**2-X**2-(1-X), X, modulus=3).is_zero
    print("Verified the paired quartic equations for arbitrary unit residue X other than 1.")


def check_walsh_descent():
    """Compare the two expansions as integer linear forms in all 16 multipliers."""
    beta = (0, 3, 2, 1)
    for word in product(range(4), repeat=4):
        a, b, c, d = word
        direct = [0]*16
        for low in product(range(2), repeat=4):
            old = tuple((x & 2) + t for x, t in zip(word, low, strict=True))
            i, j, k, l = old
            if i ^ k != j ^ l:
                continue
            walsh = sum((x & 1)*t for x, t in zip(word, low, strict=True))
            for sector in range(4):
                exponent = walsh + (sector & k).bit_count() + (beta[sector] & l).bit_count()
                direct[(i ^ k) + 4*sector] += (-1)**exponent
        folded = [0]*16
        h = a ^ c
        if h == b ^ d:
            for k, j in product(range(2), range(h & 1, 4, 2)):
                exponent = ((a ^ b) & 1)*k + (j >> 1)*(c >> 1) + (beta[j] >> 1)*(d >> 1)
                index = (h & 2)+k+4*j
                folded[index] = 4*(-1)**exponent
                assert (index & 2) == ((a ^ c) & 2)
                assert (index >> 2) % 2 == ((a ^ c) & 1)
        assert direct == folded
    print("Compared 256 Walsh tensor entries and both descent characters.")


def encode(value):
    return sum(int(c)*3**j for j, c in enumerate(value.to_list()))


def quadratic_grid():
    field = fq_default_ctx(3, 24, "a")
    ring = fq_default_poly_ctx(field)
    t = ring.gen()
    subfield = sorted((value for value, _ in (t**9-t).roots()), key=encode)
    assert len(subfield) == 9
    count = 0
    for x in subfield:
        if x**3 == x or encode(x) > encode(x**3):
            continue
        c = x.sqrt()
        for y in subfield:
            if y.is_zero():
                continue
            roots = ((t-t**3)**2-y/x**3).roots()
            assert len(roots) == 6
            v = roots[0][0]
            k, h = 1/v, c*(1/v-v)
            a, b = ((k*k-1+sign*2*h).sqrt() for sign in (1, -1))
            n = (a+b)/2
            p, other_p, q, other_q = ((2*(k+delta+sign*h)).sqrt()
                                      for delta, sign in ((1, 1), (1, -1), (-1, 1), (-1, -1)))
            coordinates = tuple(map(encode, (k, h, n, p, other_p, q, other_q)))
            verify_unit_point(field, x, y, coordinates, ((0, 1, 0), (0, 0, 1)))
            count += 1
            print("Verified the quadratic residue pair", encode(x), encode(y), flush=True)
    assert count == 24


if __name__ == "__main__":
    unit_cover()
    check_walsh_descent()
    quadratic_grid()
