"""Check Euler dimensions and Chern data of the elliptic fiber transform."""
import sympy as sp


def main():
    for degree, euler in [(3, 1), (4, 4)]:
        assert 4*degree - (degree + 8) == euler
    # Basis (rank, h, F, rank + ch2); h^2=4, hF=8, F^2=0.
    unit = sp.Matrix([1, 0, 0, 1])
    hyperplane = sp.Matrix([1, 1, 0, 3])
    input_sheaf = sp.Matrix([0, 0, 1, -4])
    universal_extension = input_sheaf + 4*unit
    cokernel = 4*hyperplane - universal_extension
    assert cokernel == sp.Matrix([0, 4, -1, 12])
    assert 4*cokernel[1]**2 + 16*cokernel[1]*cokernel[2] == 0
    assert 4*12 - (4 + 8)*4 == 0
    print('PASS: multiplication dimensions and both fiber-transform vectors')


if __name__ == '__main__':
    main()
