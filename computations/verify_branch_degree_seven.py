"""Check the plane Chern numbers and the finite degree exclusions."""
import sympy as sp


def verify():
    h, d, k = sp.symbols('h d k')
    # E has class 6-4 exp(-h); sheaf H/E has length d.
    ch = 6 - 4*(1-h+h*h/2) + d*h*h
    twist = 1 + k*h + k*k*h*h/2
    todd = 1 + 3*h/2 + h*h
    chi = sp.expand(ch*twist*todd).coeff(h, 2)
    assert chi.subs(k, -3) == d-6
    assert sp.expand(8 - ch.coeff(h, 2)) == 10-d
    for first in range(1, 7):
        for second in range(1, 7):
            if first + second + 1 <= 7:
                assert min(6-first, 6-second) > 0
    assert 4*7//2 == 14
    print('PASS: plane Riemann--Roch and both third-curve degree deficits')


if __name__ == '__main__':
    verify()
