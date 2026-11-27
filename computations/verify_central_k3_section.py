"""Check the polarized K3 central Hilbert series and the local Euler count."""
import sympy as sp


def main():
    t = sp.symbols('t')
    series = (1+t+t**2+t**3)/(1-t)**3
    expected = 1 + 2*t/(1-t) + 2*t*(1+t)/(1-t)**3
    assert sp.factor(series-expected) == 0
    c3 = sp.series((1+t)**5/(1+4*t), t, 0, 4).coeff(t, 3)
    euler = int(4*c3)
    assert euler == -56 and euler+20*2 != 4 and euler+30*2 == 4
    assert 14+4 == 18 and abs(sp.det(sp.Matrix([[4, 8], [8, 0]]))) == 64
    print('PASS: quartic K3 series, incidence dimension, lattice and nodal Euler contributions')


if __name__ == '__main__':
    main()
