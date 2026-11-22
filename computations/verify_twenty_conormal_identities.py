"""Independence of the twenty moving-point Schur identities at two rational windows."""

from itertools import product

from flint import nmod_mat

from verify_klein_elliptic_completion import rational_example
from verify_klein_full_rank_configurations import partner
from verify_klein_two_doublets import modular_points
from verify_unrestricted_doublet_deformations import K, exchange, symmetric_example


def check(w, prime):
    second_cut = exchange(w)
    first = modular_points(w, prime)
    second = partner(first, w, False, prime)
    third = partner(second, second_cut, True, prime)
    fourth = partner(third, w.T, False, prime)
    assert partner(fourth, second_cut.T, True, prime) == first

    # The row and column pairs of the second cut are (k,j) and (i,l).
    # Check both rank-three partial evaluations, not just their common zeros.
    def reduced(matrix):
        return nmod_mat(
            matrix.rows,
            matrix.cols,
            [int(x.p) * pow(int(x.q), -1, prime) % prime for x in matrix],
            prime,
        )

    matrix = reduced(second_cut)
    rows = []
    for a, b, c, d in zip(first, second, third, fourth):
        insert_b = nmod_mat(
            [[int(k == r) * b[j] for r in K] for k, j in product(K, repeat=2)], prime
        )
        insert_a = nmod_mat(
            [[int(l == r) * a[i] for r in K] for i, l in product(K, repeat=2)], prime
        )
        left = insert_b.transpose() * matrix
        right = matrix * insert_a
        assert left.rank() == right.rank() == 3
        assert (left * insert_a).rank() == 0
        assert (nmod_mat([list(c)], prime) * left).rank() == 0
        assert (right * nmod_mat([list(d)], prime).transpose()).rank() == 0
        rows.append(
            [c[k] * b[j] * a[i] * d[l] % prime for k, j, i, l in product(K, repeat=4)]
        )
    assert nmod_mat(rows, prime).rank() == len(rows) == 20
    print(
        f"Prime {prime}: twenty independent conormal products; both evaluations rank three."
    )


if __name__ == "__main__":
    check(symmetric_example(), 12619)
    check(rational_example()[0], 2657)
