"""One-return Koszul differentials from four rational primal relation windows."""
from functools import lru_cache
from itertools import product

from frobenius_tensor import nullspace_sparse, orthogonal_windows
from verify_quartic_obstruction import insertion_columns


def word_frame(windows, start, length):
    """A basis of the dual of a word quotient and its free coordinate words."""
    vectors = nullspace_sparse(insertion_columns(windows, start, length), 4**length)
    # In nullspace_sparse each free word occurs in exactly one vector with
    # coefficient 1. Other unit rows work just as well for taking coordinates.
    occurrences = {}
    for index, vector in enumerate(vectors):
        for word, value in vector.items():
            if word in occurrences or value != 1:
                occurrences[word] = None
            else:
                occurrences[word] = index
    coordinates = {index: word for word, index in occurrences.items() if index is not None}
    assert len(coordinates) == len(vectors)
    return vectors, [coordinates[index] for index in range(len(vectors))]


class CutFrames:
    def __init__(self, windows):
        self.windows = {"koszul": windows, "middle": orthogonal_windows(windows)}

    @lru_cache(maxsize=256)
    def frame(self, kind, start, length):
        return word_frame(self.windows[kind], start % 4, length)

    @lru_cache(maxsize=256)
    def split(self, kind, start, length, side):
        """Components indexed by basis vector, removed letter, shorter basis vector."""
        vectors, _ = self.frame(kind, start, length)
        _, words = self.frame(kind, (start + (side == "first")) % 4, length - 1)
        power = 4**(length - 1)
        slices = []
        for vector in vectors:
            by_letter = []
            for letter in range(4):
                entries = {}
                for index, word in enumerate(words):
                    full = letter * power + word if side == "first" else 4 * word + letter
                    if full in vector:
                        entries[index] = vector[full]
                by_letter.append(entries)
            slices.append(by_letter)
        return slices

    def complex(self, origin, i, j):
        """origin is the original vertex at position zero immediately after the cut."""
        if not all(0 <= v < 4 for v in (origin, i, j)):
            raise ValueError("cut origin and endpoints must be four-window positions")
        bases = {}
        for p, q in product(range(4-i), range(j+1)):
            r = i-j+p+q
            if r < 0:
                continue
            terms = (self.frame("koszul", origin+i, p),
                     self.frame("middle", origin+j-q, r),
                     self.frame("koszul", origin+j-q, q))
            basis = bases.setdefault(p+q, [])
            basis.extend((p, q, a, m, b) for a, m, b in
                         product(*(range(len(vectors)) for vectors, _ in terms)))
        indices = {h: {v: n for n, v in enumerate(base)} for h, base in bases.items()}
        differentials = {}
        for h in range(min(bases)+1, max(bases)+1):
            columns = []
            r = i-j+h
            for p, q, a, m, b in bases[h]:
                start = origin+j-q
                column = {}
                if p:
                    left = self.split("koszul", origin+i, p, "last")[a]
                    middle = self.split("middle", start, r, "last")[m]
                    for pieces, contractions in zip(left, middle):
                        for u, x in pieces.items():
                            for v, y in contractions.items():
                                row = indices[h-1][p-1, q, u, v, b]
                                column[row] = column.get(row, 0) + x*y
                if q:
                    right = self.split("koszul", start, q, "first")[b]
                    middle = self.split("middle", start, r, "first")[m]
                    for pieces, contractions in zip(right, middle):
                        for u, x in pieces.items():
                            for v, y in contractions.items():
                                row = indices[h-1][p, q-1, a, v, u]
                                column[row] = column.get(row, 0) + (-1)**p*x*y
                columns.append({row: value for row, value in column.items() if value})
            differentials[h] = columns
        return {h: len(base) for h, base in bases.items()}, differentials
