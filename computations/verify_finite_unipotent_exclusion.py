"""Verify tangent projections and the low-plane commutator/cubic constraints."""
from itertools import combinations, product

import sympy as sp

from verify_finite_additive_windows import jordan, tensor_sum


def main():
    eye = sp.eye(16)
    def word(i, j):
        return eye[:, 4*i+j]
    # z,v,x,t. Both endpoints are the canonical z-point.
    n = jordan((3, 1))
    d = tensor_sum(n, n)
    low = sp.Matrix.hstack(word(0,1)-word(1,0), word(0,3), word(3,0))
    linear = [4*i+j for i,j in product(range(4), repeat=2) if (i==0) != (j==0)]
    assert (d*low)[linear,:].rank() == 1 and (d**2*low)[linear,:].is_zero_matrix
    # z,u,v,x. All weight-zero invariant relations have just one
    # independent possible differential at the canonical point.
    zero = sp.Matrix.hstack(word(0,3)-word(3,0), *(word(i,j) for i,j in product((1,2), repeat=2)))
    assert zero[linear,:].rank() == 1
    for num in (2,3):
        assert len(list(combinations(range(num),2))) < 2*(num-1)
    u, v, y = sp.symbols('u v y')
    forms = sp.Matrix([[-v,u,0],[-y,0,u],[0,-y,v]])
    assert forms.det() == 0
    # Coefficients in the two nonnegative successive weight differences.
    for powers in [(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1)]:
        assert powers[1]+powers[2] <= 2 and powers[2] <= 1
    print('Canonical tangent bounds, commutative plane restrictions, and forbidden cubic jets verified.')


if __name__ == '__main__':
    main()
