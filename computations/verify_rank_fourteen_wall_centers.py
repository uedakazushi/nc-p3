"""Check center base change on every ordinary rank-fourteen geometric residue at 101."""

from concurrent.futures import ProcessPoolExecutor
from functools import cache

from flint import nmod_poly

from rank_fourteen_function_field import one
from verify_klein_center_lift import add_multiple, central_equations, field_paths
from verify_rank_fourteen_geometric_reduction import rational_windows, test_exception


def verify():
    bad = set()

    @cache
    def record(coefficients):
        polynomial = nmod_poly(list(coefficients), 101)
        assert polynomial
        bad.update(tuple(map(int, f)) for f, _ in polynomial.factor()[1])

    def record_reduction(polynomial):
        record(tuple(int(c) % 101 for c in polynomial))

    def reduce(rows):
        # The product, up to sign, is the actual minor of the original rows and
        # pivot columns. Cramer's rule extends the complete chart off its zero set.
        pivots = {}
        determinant = one
        for original in rows:
            row = original.copy()
            for j in sorted(pivots):
                if j in row:
                    add_multiple(row, pivots[j], -row[j])
            if row:
                j = min(row)
                determinant *= row[j]
                inverse = row[j].inverse()
                pivots[j] = {k: v*inverse for k, v in row.items()}
        for j in sorted(pivots, reverse=True):
            for k in list(pivots[j]):
                if k > j and k in pivots:
                    add_multiple(pivots[j], pivots[k], -pivots[j][k])
        record_reduction(determinant.norm_numerator())
        record_reduction(determinant.denominator)
        return pivots

    windows = rational_windows()
    for side in windows:
        for _, row in side:
            for c in row.values():
                record_reduction(c.denominator)
    basis, image = field_paths(windows, one, reduce)
    for p in range(4):
        assert [sum(len(basis[p, n, h]) for h in range(4)) for n in range(6)] == [
            1, 4, 10, 20, 35, 56
        ]
    for h, expected in enumerate((3, 2, 0, 0)):
        rows, width = central_equations(basis, image, h)
        assert width - len(reduce(rows)) == expected
    assert (len(bad), max(len(f)-1 for f in bad)) == (33, 12)
    with ProcessPoolExecutor(max_workers=4) as pool:
        cases = [(f, i) for f in sorted(bad) for i in (10, 91)]
        ordinary = sum(pool.map(test_exception, cases))
    print("33 residual factors, 2 embeddings;", ordinary, "ordinary cases certified.")


if __name__ == "__main__":
    verify()
