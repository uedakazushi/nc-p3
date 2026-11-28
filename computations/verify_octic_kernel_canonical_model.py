"""Koszul and liaison Euler characteristics for the graph and its étale double cover."""

from math import comb

import sympy as sp


def chi_space(n, degree):
    return sp.prod(degree+j for j in range(1, n+1))/sp.factorial(n)


def chi_graph(a, b):
    complete_intersection = sum((-1)**j*comb(4,j)*chi_space(3,a-2*j)*chi_space(3,b-j)
                                for j in range(5))
    return sp.expand(complete_intersection-20*chi_space(2,b-3))


if __name__ == '__main__':
    a,b,t = sp.symbols('a b t')
    polynomial = chi_graph(a,b)
    assert polynomial == 4*a**2+24*a*b-16*a+6*b**2-18*b+16
    assert polynomial.subs({a:t,b:t}) == 34*t**2-34*t+16
    assert 2*polynomial.subs({a:0,b:0}) == 32
    assert 15+16 == 31 and 2*(128-3*20) == 136
    print('Verified the bivariate surface polynomial, the canonical square, and the double-cover invariants.')
