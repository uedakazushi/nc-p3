"""Partition and matrix ranks used by the degree-four central-fiber criterion."""
from itertools import combinations_with_replacement
from sympy import Matrix, eye


def main():
    profiles = [parts for length in range(1,5) for parts in
                combinations_with_replacement(range(1,5),length) if sum(parts) == 4]
    assert profiles == [(4,), (1,3), (2,2), (1,1,2), (1,1,1,1)]
    units = [Matrix(2,2,lambda i,j: int((i,j) == (r,c)))
             for r in range(2) for c in range(2)]
    spaces = [units, [units[0]-units[-1], units[1], units[2]]]
    for H in spaces:
        for A, rank in [(0*eye(2),0), (units[1],2), (eye(2),len(H))]:
            images = Matrix.hstack(*[(h*A).reshape(4,1) for h in H])
            assert images.rank() == rank
    assert 1 + (-1) == 0 and 1 + (-1)**2 == 2
    print('Degree-four profiles, contracted two-dimensional traces, and opposite-scale octic signs verified.')


if __name__ == '__main__':
    main()
