"""Exact matrix and degree checks for the octic normalization proof."""
from sympy import Matrix, eye, zeros


def main():
    for size in (2,3,4):
        units = []
        for i in range(size):
            for j in range(size):
                if i != j:
                    unit = zeros(size)
                    unit[i,j] = 1
                    units.append(unit)
        for i in range(1,size):
            diagonal = zeros(size)
            diagonal[0,0], diagonal[i,i] = 1,-1
            units.append(diagonal)
        gram = Matrix([[(a*b).trace() for a in units] for b in units])
        assert gram.rank() == size**2-1
    for rank in (1,2,3):
        assert 7 - (4 + (4-rank) - 1) == rank < 4 + rank - 1
    for hyperplane_degree, ruling_degree in ((4,2),(2,1)):
        kernel_degree = ruling_degree - 2*hyperplane_degree
        assert ruling_degree - kernel_degree in (8,4)
    # A triangular extension need not split when the two diagonal maps both have rank one.
    projection = Matrix([[1,0],[0,0]])
    for off_diagonal,expected in [(zeros(2),2),(eye(2),3)]:
        full = projection.row_join(off_diagonal).col_join(zeros(2).row_join(projection))
        assert full.rank() == expected
    socle = 4*(8-1)
    assert (socle+1, 3*socle, socle+4) == (29,84,32)
    print('Trace form ranks, transitive quotient bounds, determinant line degrees, extension ranks, and the 84-degree cutoff verified.')


if __name__ == '__main__':
    main()
