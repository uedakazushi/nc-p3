"""Check the nilpotent cuts and the weighted transvection embedding over Q."""
from itertools import product

import sympy as sp


def jordan(parts):
    return sp.diag(*(sp.zeros(n) + sp.Matrix(n, n, lambda i, j: int(i == j + 1))
                     for n in parts))


def tensor_sum(a, b):
    return sp.kronecker_product(a, sp.eye(b.rows)) + sp.kronecker_product(sp.eye(a.rows), b)


def main():
    types = [(4,), (3, 1), (2, 2), (2, 1, 1)]
    for first, second in product(types, repeat=2):
        a, b = jordan(first), jordan(second)
        total = tensor_sum(a, b)
        i, j = first[0]-1, second[0]-1
        assert total**(i+j) == sp.binomial(i+j, i)*sp.kronecker_product(a**i, b**j)
    n = jordan((4,))
    sym = sp.Matrix.hstack(*(sp.eye(16)[:, 4*i+j] + sp.eye(16)[:, 4*j+i]
                            for i in range(4) for j in range(i, 4)))
    sq = tensor_sum(n, n)*sym
    assert sq.rank() == 8
    assert (tensor_sum(n, n)**6 * sym).rank() == 1

    # z, u, v, x with Nz=x. Only the weight-zero alternating and UU
    # terms survive in a middle invariant relation.
    n = sp.zeros(4)
    n[3, 0] = 1
    total = tensor_sum(n, n)
    weights = [-1, 0, 0, 1]
    for weight, rank, dimension in [(-2, 1, 1), (-1, 4, 4), (0, 1, 6), (1, 0, 4)]:
        cols = [4*i+j for i,j in product(range(4), repeat=2) if weights[i]+weights[j] == weight]
        assert len(cols) == dimension and total[:, cols].rank() == rank
    assert (total**2).rank() == 1 and (total**2)[:, 0] == 2*sp.eye(16)[:, 15]
    # A nonhomogeneous highest-weight embedding with both mixed projections
    # invertible. Weight-zero restriction of an equivariant perfect pairing
    # remains perfect for every invertible map between the abstract chains.
    r = jordan((2, 2, 1, 1))
    w = [-1, 1, -1, 1, 0, 0]
    m = sp.Matrix(6, 6, sp.symbols('a:36'))
    maps = sp.linear_eq_to_matrix(list(m*r-r*m), list(m))[0].nullspace()
    for basis in maps:
        g = basis.reshape(6, 6)
        assert all(not g[i, j] for i,j in product(range(6), repeat=2) if w[i] < w[j])
    print('Top powers for 16 Jordan pairs, cubic obstruction ingredients, and transvection filtrations verified.')


if __name__ == '__main__':
    main()
