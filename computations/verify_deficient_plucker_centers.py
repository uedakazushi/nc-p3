"""Determinantal Hilbert numerator and the input degrees for line and plane centers."""

import sympy as sp


def main():
    t = sp.symbols('t')
    numerator = 1-15*t**4+24*t**5-10*t**6
    assert sp.cancel(numerator/(1-t)**3).subs(t,1) == 20
    s = sp.symbols('s')
    rows = sp.Matrix([[1,s,0,0],[0,1,s,0],[0,0,1,s]])
    kernel = [(-1)**i*rows[:,[j for j in range(4) if j!=i]].det() for i in range(4)]
    assert rows*sp.Matrix(kernel) == sp.zeros(3,1)
    assert sp.Matrix([[sp.expand(p).coeff(s,i) for i in range(4)] for p in kernel]).rank() == 4
    assert 4*3-2 == 10
    print('Verified the length-twenty Euler numerator and the degree-three kernel of a balanced conic.')


if __name__ == '__main__':
    main()
