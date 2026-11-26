"""Check the stress formulas for punctured secant configurations over Q."""
from itertools import combinations
from sympy import Matrix, zeros


def main():
    edges = list(combinations(range(6), 2))
    a = Matrix([1, 2, 3, 5, 7, -18])
    z = Matrix.eye(5).row_join(-Matrix.ones(5, 1))
    u = Matrix.hstack(*(z.col(i)/a[i] for i in range(6)))
    for holes in ([], [(0, 1)], [(0, 1), (0, 2)], [(0, 1), (2, 3)]):
        used = [ij for ij in edges if ij not in holes]
        forms = Matrix([[int(k in ij) for k in range(6)] for ij in edges])
        conditions = [a.T] + [forms.row(edges.index(ij)) for ij in holes]
        params = Matrix.hstack(*Matrix.vstack(*conditions).nullspace())
        coeffs = Matrix.vstack(*(a[i]*a[j]*forms.row(edges.index((i, j)))
                                 for i, j in used))*params
        squares = []
        for i, j in used:
            d = u.col(i)-u.col(j)
            squares.append((d*d.T).reshape(25, 1))
        evals = Matrix.hstack(*squares)
        assert evals*coeffs == zeros(25, params.cols)
        assert coeffs.rank() == len(evals.nullspace()) == len(used)-10
    print("PASS: complete square dependencies of 13, 14, 15 labeled secant intersections")


if __name__ == "__main__":
    main()
