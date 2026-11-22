"""Recover the two-dimensional letter blocks at all four edges of two full-rank windows."""

import sympy as sp
from flint import nmod_mat

from verify_companion_rank_two_sections import PRIME, orbit_types, section
from verify_klein_elliptic_completion import gram_twist


def check():
    q = sp.Rational
    samples = [
        (
            sp.Matrix([[q(36, 17), 0], [9, 4]]),
            sp.Matrix([[q(36, 19), q(36, 37)], [q(36, 5), -12]]),
            3,
            2,
            [(1, 2, 3, 4, 5, 15), (1, 2, 3, 4), (1, 2, 3, 4, 5), (1, 2, 4, 5, 6, 10)],
        ),
        (
            sp.Matrix([[1, 1], [1, 2]]),
            sp.Matrix([[-q(13, 14), q(1, 13)], [-14, -q(16, 11)]]),
            -2,
            -2,
            [(1, 3, 4, 5), (1, 2, 3, 4), (1, 2, 3, 5, 6), (1, 2, 4, 5, 7, 12)],
        ),
    ]
    for A, B, z, l, seeds in samples:
        w, opp, rev = gram_twist(A, B, sp.Integer(z), sp.Integer(l))
        assert opp.rank() == 16
        for pos, (mat, values) in enumerate(zip([w, rev, w.T, rev.T], seeds)):
            columns = mat.columnspace()
            rels = [
                [
                    int(sp.numer(c)) * pow(int(sp.denom(c)), -1, PRIME) % PRIME
                    for c in vec
                ]
                for vec in columns
            ]
            assert nmod_mat(rels, PRIME).rank() == 6
            types = set.intersection(
                *(orbit_types(section(rels, seed)) for seed in values)
            )
            assert types == {((1, 4), (1, 16)), ((1, 20),)}
            thick = [
                v.reshape(4, 4)
                for v in columns
                if not any(v[4 * i + j] for i in range(2) for j in range(2, 4))
            ]
            assert len(thick) == 4
            system = sp.Matrix.vstack(
                *[
                    sp.kronecker_product(sp.eye(4), m.T).row_join(
                        -sp.kronecker_product(m, sp.eye(4))
                    )
                    for m in thick
                ]
            )
            assert system.rank() == 30
            print(
                f"z={z}, slope={l}, position={pos}: 4+16, two-dimensional common endomorphism algebra."
            )


if __name__ == "__main__":
    check()
