"""Closed points of degrees two through four on the negative unit residue line."""

from itertools import product

import sympy as sp
from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx

from verify_rank_fifteen_unit_residues import (
    Scalar, graded_profile, pattern, tensor_relations, verify_unit_point,
)


def negative_x_coordinates():
    w, g, k = sp.symbols("w g k")
    s = sp.I+w
    n = g*(1-s**2)/(2*sp.I+w)
    other = 2*g*s/(2*sp.I+w)
    h = n*other
    assert sp.cancel(n**2+other**2-(g*w)**2) == 0
    assert (n.subs(w, 0), other.subs(w, 0)) == (-sp.I*g, g)
    for sign in (1, -1):
        p_squared = 2*(k+1+sign*h)
        q_scaled_squared = 2*w**2/(k+1)+sign*2*h/g**2
        assert sp.simplify(p_squared.subs({k: 1, w: 0})-(4-sign*2*sp.I*g**2)) == 0
        assert sp.simplify(q_scaled_squared.subs({k: 1, w: 0})+sign*2*sp.I) == 0
    t = sp.symbols("t")
    a, b = (1-t**2)/(1+t**2), 2*t/(1+t**2)
    g_squared = (2*(1-sp.I))**2*(t-t**3)/(1+t**2)**2
    for sign in (1, -1):
        assert sp.cancel((2*(a-sign*b))**2-(4-sign*2*sp.I*g_squared)) == 0
    assert sp.cancel(-g_squared**2/16-4*(t-t**3)**2/(1+t**2)**4) == 0
    assert sp.cancel(b/(1+a)-t) == 0

    field = fq_default_ctx(3, 2, "a")
    ring = fq_default_poly_ctx(field)
    t = ring.gen()
    i = (-field.one()).sqrt()
    plus, minus = 4+2*i*t**2, 4-2*i*t**2
    assert plus.gcd(minus).degree() == 0
    assert all(f.gcd(f.derivative()).degree() == 0 for f in (plus, minus))
    assert plus(i*t) == minus and minus(i*t) == plus
    assert not (t**4+4) % plus and not (t**4+4) % minus


def negative_x_prime_slopes():
    field = fq_default_ctx(3, 8, "a")
    ring = fq_default_poly_ctx(field)
    i = (-field.one()).sqrt()
    for slope in (1, -1):
        g = (ring.gen()**4+16*slope).roots()[0][0]
        n, h = -i*g, -i*g*g
        p, u = (4+2*h).sqrt(), (4-2*h).sqrt()
        q, v = (1-i)*g, (1+i)*g
        assert g**4 == field(-16*slope) and q*q == 2*h and v*v == -2*h
        gamma = pattern(field.one(), field.zero(), n/2, g/2,
                        (u-p)/4, (u+p)/4, (v-q)/4, (v+q)/4)
        entries = tuple(Scalar(field.zero() if c == 0 else c) for c in gamma)
        windows = tensor_relations(entries, Scalar(field.zero()))
        assert graded_profile(windows, Scalar(field.one()), ((0, 1, 0), (0, 0, 1)))


def verify():
    negative_x_coordinates()
    negative_x_prime_slopes()
    prime = fmpz_mod_poly_ctx(3)
    for degree, expected in ((2, 3), (3, 8), (4, 18)):
        field = fq_default_ctx(3, 12*degree, "a")
        ring = fq_default_poly_ctx(field)
        t = ring.gen()
        i = (-field.one()).sqrt()
        count = 0
        for coeff in product(range(3), repeat=degree):
            polynomial = [*coeff, 1]
            if not prime(polynomial).is_irreducible():
                continue
            y = ring(polynomial).roots()[0][0]
            v = ((t-t**3)**2+y).roots()[0][0]
            k = 1/v
            h = i*(k-v)
            assert -(v-v**3)**2 == y
            a, b = ((k*k-1+sign*2*h).sqrt() for sign in (1, -1))
            n = (a+b)/2
            p, u, q, w = ((2*(k+delta+sign*h)).sqrt()
                       for delta, sign in ((1, 1), (1, -1), (-1, 1), (-1, -1)))
            encoded = tuple(sum(int(c)*3**j for j, c in enumerate(value.to_list()))
                            for value in (k, h, n, p, u, q, w))
            verify_unit_point(field, -field.one(), y, encoded, ((0, 1, 0), (0, 0, 1)))
            count += 1
            print("Verified negative unit residue orbit", polynomial, flush=True)
        assert count == expected


if __name__ == "__main__":
    verify()
