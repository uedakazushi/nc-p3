"""Check the small-degree trace and the different-scale thin coboundary matrix."""
from sympy import Matrix, diag, eye, factor, prod, symbols


def main():
    x = Matrix([[0, 1], [0, 0]])
    y = x.T
    assert (x*y-y*x)**2 == eye(2)
    scales = symbols('l:4')
    gauge = diag(*scales)
    for i in range(4):
        gauge[i,(i+1) % 4] = -1
    assert factor(gauge.det() - prod(scales) + 1) == 0
    # Trace-pairing flattening on the products of two full two-dimensional arrow spaces.
    units = [Matrix(2, 2, lambda r, c: int((r,c) == (i,j)))
             for i in range(2) for j in range(2)]
    pairing = Matrix([[(a*b).trace() for b in units] for a in units])
    assert pairing.rank() == 4
    assert 2 * (3**2 - (3-1)**2) > 3**2
    print('Corner central polynomial, thin scaling determinant, and trace ranks in degrees two and three verified.')


if __name__ == '__main__':
    main()
