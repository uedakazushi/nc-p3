"""Check both generic flag Fourier transforms and their slope propagation."""
import sympy as sp


def main():
    assert -1 - (4*(-1) - (-1)) == 2
    first = sp.Matrix([[1, 3, 5, 7], [0, -8, -16, -24]])
    second = sp.Matrix([[3, 1, -1, -3], [32, 8, -16, -40]])
    for degree, classes, sign, expected in [
            (-4, first, -1, [-1, -3, -5, -7]),
            (12, second, 1, [3, 1, -1, -3])]:
        transform = sp.Matrix([[0, 1], [-1, degree]])
        bundles = sign*transform.inv()*classes
        assert list(bundles.row(0)) == [4]*4
        assert list(bundles.row(1)) == expected
        for i in range(4):
            for j in range(4):
                slope_difference = (bundles[1, j] - bundles[1, i])/4
                assert slope_difference == sp.Rational(i-j, 2)
                assert slope_difference + 2 > 0
    print('PASS: determinant pencil, generic ranks, signed degrees and both slope increments')


if __name__ == '__main__':
    main()
