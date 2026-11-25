"""Exact ternary boundary eigenranks, node evaluations and mixed-coefficient lifting."""

import sympy as sp
from ternary_ore_parameters import kernel, presentation_pencil


def relations(zeta, coefficients):
    e = sp.eye(9)
    return sp.Matrix.hstack(
        *[
            coefficients[j] * e[:, 4 * j]
            + e[:, 3 * ((j + 1) % 3) + (j + 2) % 3]
            - zeta * e[:, 3 * ((j + 2) % 3) + (j + 1) % 3]
            for j in range(3)
        ]
    )


def main():
    for zeta in [-1, sp.I]:
        expected = {1: 3, -1: 6} if zeta == -1 else {-1: 3, sp.I: 3, -sp.I: 3}
        for h in range(3):
            cs = [1] * h + [0] * (3 - h)
            left, right, _, lift = presentation_pencil(relations(zeta, cs))
            spaces = []
            for mu, dim in expected.items():
                eig = kernel(left + mu * right)
                # The 3-dimensional gauge is exactly L_v - mu R_v.
                assert eig.cols == dim + 3
                spaces.append(eig)
            if h == 1:
                weights = (1, 0, 2)
                quadratic_weights = [
                    next(
                        weights[r // 3] + weights[r % 3] for r in range(9) if lift[r, a]
                    )
                    for a in range(6)
                ]
                for weight in (-1, 3):
                    cols = [
                        i
                        for i in range(18)
                        if quadratic_weights[i % 6] - weights[i // 6] == weight
                    ]
                    eigendims = tuple(
                        len(cols) - (left + mu * right)[:, cols].rank()
                        for mu in expected
                    )
                    assert eigendims == (
                        (1, 1)
                        if zeta == -1
                        else ((1, 0, 1) if weight == -1 else (1, 1, 0))
                    )
            assert left.row_join(right).rank() == 18
            all_eigen = sp.Matrix.hstack(*spaces)
            for j in range(h, 3):
                # At vertex j the two other components of d(p,p) must not both vanish.
                evaluations = sp.zeros(2, 18)
                for i, a in enumerate(b for b in range(3) if b != j):
                    evaluations[i, 6 * a : 6 * a + 6] = lift[4 * j, :]
                assert (evaluations * all_eigen).rank() == 2
            print("PASS", zeta, cs, "eigenspace and vertex ranks")
    print(
        "PASS: full boundary eigenvectors lift, including every quantum-base correction"
    )


if __name__ == "__main__":
    main()
