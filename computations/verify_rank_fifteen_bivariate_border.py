"""Thirteen-degree four-center certificates over the two-variable residue field."""

from verify_rank_fifteen_bivariate_seventh_cut import (
    clear_denominators,
    eliminate,
    fraction,
    pivot_degree,
)
from verify_rank_fifteen_bivariate_word_paths import (
    adjust,
    factors,
    field,
    record,
    tensor,
)
from verify_rank_fifteen_fixed_center_obstruction import universal_sections
from verify_rank_fifteen_center_border import border_data
from central_border_paths import central_border

EXPECTED = (4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1, 0)


def reduce_words(rows, n):
    if n > 8:
        return {}, [r for r in rows if r]
    pending = [clear_denominators(r) for r in rows if r]
    pivots = []
    residuals = []
    while pending:
        widths = [sum(len(w) == n for w in r) for r, _ in pending]
        if 0 in widths:
            residuals.extend(
                {w: fraction(z) * scale for w, z in r.items()}
                for (r, scale), width in zip(pending, widths, strict=True)
                if not width
            )
            pending = [
                entry for entry, width in zip(pending, widths, strict=True) if width
            ]
            widths = [width for width in widths if width]
            if not pending:
                break
        width = min(widths)
        _, column, index = min(
            (pivot_degree(z, scale, True), w, i)
            for i, (r, scale) in enumerate(pending)
            if widths[i] == width
            for w, z in r.items()
            if len(w) == n
        )
        pivot, scale = pending.pop(index)
        record((scale * fraction(pivot[column])).n)
        pivots.append((column, pivot))
        for i, (other, multiplier) in enumerate(pending):
            if column in other:
                pending[i] = eliminate(other, multiplier, pivot, column)
        pending = [(r, scale) for r, scale in pending if r]
    reduced = {}
    for column, row in reversed(pivots):
        for w in list(row):
            if w in reduced:
                row, _ = eliminate(row, field.one, reduced[w], w)
        reduced[column] = row
    result = {
        j: {w: fraction(z, row[j]) for w, z in row.items()}
        for j, row in reduced.items()
    }
    return result, residuals


def independent_sources(rows):
    pending = [(i, *clear_denominators(row)) for i, row in enumerate(rows) if row]
    selected = []
    while pending:
        width = min(len(r) for _, r, _ in pending)
        _, column, index = min(
            (pivot_degree(z, scale, True), w, i)
            for i, (_, r, scale) in enumerate(pending)
            if len(r) == width
            for w, z in r.items()
        )
        source, pivot, scale = pending.pop(index)
        record((scale * fraction(pivot[column])).n)
        selected.append(rows[source])
        for i, (idx, row, multiplier) in enumerate(pending):
            if column in row:
                pending[i] = (idx, *eliminate(row, multiplier, pivot, column))
        pending = [item for item in pending if item[1]]
    return selected


def certificate(centers, axis=0, start=0):
    relations, weighted = border_data(centers, tensor, adjust, start)
    slope = [field.one if j == axis else field.zero for j in range(5)]
    for n, dims in central_border(
        relations,
        weighted,
        slope,
        field.one,
        reduce_words=reduce_words,
        select=independent_sources,
    ):
        assert sum(dims) == EXPECTED[n - 1]
        print(n, dims, flush=True)
    return set(factors)


if __name__ == "__main__":
    divisors = certificate(universal_sections(range(5)))
    assert len(divisors) == 171
    print("PASS: two-variable exchanged four-cut through thirteen, 171 pivot factors")
