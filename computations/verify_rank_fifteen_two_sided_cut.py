"""Saturate two-sided path overlaps in the four-central-quartic quotient."""

from itertools import product

from verify_klein_center_lift import add_multiple
from verify_rank_fifteen_variable_zero_disk import center_data
from verify_rank_fourteen_all_ordinary import DIMENSIONS, saturated_projection


def initial_layers(jet, build_relations):
    _, basis, image, centers = center_data(15, jet, build_relations)
    precision = jet.precision
    assert precision == 11
    dims, left, right, words = {}, {}, {}, {}
    for n in range(6):
        projections = {}
        for start, h in product(range(4), repeat=2):
            jet.precision = precision
            paths = basis[start, n, h]
            rows = []
            if n >= 4:
                for g, center in centers[start]:
                    for tail in basis[start, n-4, g^h]:
                        row = {}
                        for w, c in center.items():
                            add_multiple(row, image(start, w+tail), c)
                        rows.append(row)
            free, pi, loss = saturated_projection(rows, len(paths), DIMENSIONS[n][h], jet)
            assert loss == (4 if n == 5 else 0)
            words[start, n, h] = [paths[j] for j in free]
            projections[start, h] = pi
            dims[start, n, h] = len(free)
        if n == 0:
            continue
        jet.precision = precision-(4 if n == 5 else 0)
        for start, h in product(range(4), repeat=2):
            def coords(w):
                result = {}
                for j, coefficient in image(start, w).items():
                    add_multiple(result, projections[start, h][j], coefficient)
                return result

            for a in range(4):
                right[start, n, h, a] = [coords(w+(a,)) for w in words[start, n-1, h^a]]
                left[start, n, h, a] = [coords((a,)+w) for w in words[(start+1)%4, n-1, h^a]]
    assert jet.precision == 7
    return dims, left, right


def overlap_rows(n, start, h, dims, left, right):
    following = (start+1) % 4
    ends, beginnings = [], []
    width = 0
    for p, offsets in ((start, ends), (following, beginnings)):
        for a in range(4):
            offsets.append(width)
            width += dims[p, n-1, h^a]
    rows = []
    for a, b in product(range(4), repeat=2):
        for j in range(dims[following, n-2, h^a^b]):
            row = {ends[b]+k: v for k, v in left[start, n-1, h^b, a][j].items()}
            row.update({beginnings[a]+k: -v for k, v in right[following, n-1, h^a, b][j].items()})
            rows.append(row)
    return ends, beginnings, width, rows


def store_projection(n, start, h, dims, left, right, ends, beginnings, free, projection):
    following = (start+1) % 4
    dims[start, n, h] = len(free)
    for a in range(4):
        right[start, n, h, a] = projection[ends[a]:ends[a]+dims[start, n-1, h^a]]
        left[start, n, h, a] = projection[beginnings[a]:beginnings[a]+dims[following, n-1, h^a]]


def overlap(n, start, h, dims, left, right, jet):
    ends, beginnings, width, rows = overlap_rows(n, start, h, dims, left, right)
    free, projection, loss = saturated_projection(rows, width, DIMENSIONS[n][h], jet)
    store_projection(n, start, h, dims, left, right, ends, beginnings, free, projection)
    return loss


def verify(jet, build_relations):
    dims, left, right = initial_layers(jet, build_relations)
    for n in range(6, 14):
        losses = [overlap(n, p, h, dims, left, right, jet)
                  for p, h in product(range(4), repeat=2)]
        assert losses == [2 if n in (7, 9) or (n == 6 and h >= 2) else 0
                          for p, h in product(range(4), repeat=2)]
        jet.precision -= max(losses)
        assert jet.precision > 0
        print("Length", n, "central quotient bound", DIMENSIONS[n],
              "remaining precision", jet.precision, flush=True)
