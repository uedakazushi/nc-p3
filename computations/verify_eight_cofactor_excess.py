"""Numerical identities for the excess-curve obstruction on the compactified elliptic pencil."""

from itertools import product

import sympy as sp


def main():
    h=sp.symbols('h')
    # The Chern polynomial of the quotient of four O(H)'s by a numerically trivial line.
    z=sp.symbols('z')
    virtual_length = sp.expand((1+z*h)**4).coeff(z,3).subs(h**3,4)
    assert virtual_length == 16
    degree_A_dual = 4
    degree_sym_square = 4*degree_A_dual
    assert 3*degree_sym_square-4*2*degree_A_dual == 16
    for j,extra in product((1,2),range(4)):
        degree = j+extra
        min_chi = 2-j
        bound = int(virtual_length)+2*degree-4*j-min_chi
        assert bound <= (19 if j==1 else 18)
    print('Verified the two length-sixteen formulas and the cubic bounds nineteen and eighteen.')


if __name__ == '__main__':
    main()
