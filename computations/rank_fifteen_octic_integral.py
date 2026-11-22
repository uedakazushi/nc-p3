"""Integral equations of a genus-three four-cut in characteristic zero."""

from rank_fifteen_octic_field import UNIT, windows
from verify_rank_fifteen_octic_center import chosen_rows
from verify_klein_center_lift import add_multiple, field_paths, central_equations
from verify_rank_fifteen_octic_boundary import clear_coordinates
from functools import cache


def projector(rows, central=False):
    """Compute the exact column projection with the verified split pivot choices."""
    pivots = {}
    for row, k in chosen_rows(rows, central):
        for a, prev in pivots.items():
            if a in row:
                add_multiple(row, prev, -row.pop(a))
        c = row.pop(k).inverse()
        pivots[k] = {a: v * c for a, v in row.items()}
    for k in reversed(pivots):
        row = pivots[k]
        for a in list(row):
            if a in pivots:
                add_multiple(row, pivots[a], -row.pop(a))
    for k, row in pivots.items():
        row[k] = UNIT
    return pivots


@cache
def integral_presentation(color):
    """Clear each central form by one scalar common to all four positions."""
    w = windows(color)
    bs, img = field_paths(w, UNIT, projector)
    assert all(
        sum(len(bs[p, n, h]) for h in range(4)) == d
        for p in range(4)
        for n, d in enumerate((1, 4, 10, 20, 35, 56))
    )
    cs = [[] for p in range(4)]
    for h, sel in ((0, (0, 2)), (1, (0, 1))):
        rows, size = central_equations(bs, img, h)
        pivs = projector(rows, True)
        free = [k for k in range(size) if k not in pivs]
        assert len(free) == 3 - h
        for j in sel:
            col = free[j]
            vec = {col: UNIT}
            vec.update({k: -row[col] for k, row in pivs.items() if col in row})
            parts = {}
            off = 0
            for p in range(4):
                words = bs[p, 4, h]
                parts.update(
                    {
                        (p, *w): vec[off + i]
                        for i, w in enumerate(words)
                        if off + i in vec
                    }
                )
                off += len(words)
            assert off == size
            polys = clear_coordinates(parts)
            for p in range(4):
                cs[p].append((h, {k[1:]: v for k, v in polys.items() if k[0] == p}))
    rs = [[(h, clear_coordinates(row)) for h, row in side] for side in w]
    return rs, cs
