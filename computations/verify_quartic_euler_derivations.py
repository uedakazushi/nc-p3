"""Upper bounds for normalized degree-four derivations at the full-rank AS sample."""

from itertools import product

from flint import nmod_mat

from verify_klein_elliptic_completion import four_relations, rational_example
from verify_klein_period_four_center import K, paths


def check():
    prime = 1009
    w, _, rev = rational_example()
    relations = four_relations(w, rev)
    basis, image = paths(relations, 6)
    assert all(
        sum(len(basis[p, n, h]) for h in K) == size
        for p in K
        for n, size in [(4, 35), (5, 56), (6, 84)]
    )

    def residue(value):
        return int(value.p) * pow(int(value.q), -1, prime) % prime

    for h in K:
        # A derivation vanishing on vertex idempotents is determined on arrows.
        variables = [
            (p, i, word) for p, i in product(K, repeat=2) for word in basis[p, 5, i ^ h]
        ]
        columns = {v: j for j, v in enumerate(variables)}
        equations = []
        for p in K:
            for grade, relation in relations[p]:
                block = [[0] * len(columns) for _ in basis[p, 6, grade ^ h]]
                for ij, value in enumerate(relation):
                    if not value:
                        continue
                    i, j = divmod(ij, 4)
                    c = residue(value)
                    for place, arrow, suffix in [(p, i, True), ((p + 1) % 4, j, False)]:
                        for word in basis[place, 5, arrow ^ h]:
                            col = columns[place, arrow, word]
                            term = word + (j,) if suffix else (i,) + word
                            for row, coeff in zip(block, image(p, term)):
                                row[col] = (row[col] + c * residue(coeff)) % prime
                equations.extend(block)
        rank = nmod_mat(equations, prime).rank()
        quartics = sum(len(basis[p, 4, h]) for p in K)
        central = (3, 2, 0, 0)[h]
        assert len(variables) - rank == quartics == (36, 40, 32, 32)[h]
        print(
            f"Character {h}: {len(variables)} arrow coefficients, rank {rank}, "
            f"{quartics - central} inner and {central} Euler directions."
        )


if __name__ == "__main__":
    check()
