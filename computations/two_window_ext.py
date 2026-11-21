"""Four-term Ext complex for a right module on the forward-labeled cut B."""
from fractions import Fraction
from functools import lru_cache
from itertools import product

from frobenius_tensor import orthogonal_windows, reconstructed_windows
from verify_quartic_obstruction import relation_matrix
from general_cut_tor import word_frame


def pbw_overlap_counts(windows):
    """Leading-word overlap counts in the supplied ordered arrow bases (left lex)."""
    leading = []
    for relations in windows:
        if len(relations) != 6:
            raise ValueError("expected six quadratic relations per window")
        _, pivots = relation_matrix(relations).T[:, ::-1].rref()
        if len(pivots) != 6:
            raise ValueError("dependent quadratic relations")
        leading.append({divmod(15-pivot, 4) for pivot in pivots})
    return [sum((a, b) in leading[i] and (b, c) in leading[(i+1) % 4]
                for a, b, c in product(range(4), repeat=3)) for i in range(4)]


def complete_windows(first, second):
    """Use the cubic overlap coevaluation, preserving the two original windows."""
    if any(len(window) != 6 or relation_matrix(window).rank() != 6 for window in (first, second)):
        raise ValueError("two independent six-dimensional relation windows are required")
    dual = orthogonal_windows([first, second])
    overlap, _ = word_frame(dual, 0, 3)
    if len(overlap) != 4:
        raise ValueError("the cubic overlap must have dimension four")
    tensor = {4*word + i: value for i, vector in enumerate(overlap)
              for word, value in vector.items()}
    returning = reconstructed_windows(tensor)
    return [first, second, returning[2], returning[3]]


def add_scaled(vector, other, scale):
    for index, value in other.items():
        updated = vector.get(index, 0) + scale*value
        if updated:
            vector[index] = updated
        else:
            vector.pop(index, None)


def quotient_images(size, relations):
    """Sparse quotient map on every standard basis vector."""
    pivots = {}
    for relation in relations:
        row = dict(relation)
        while row:
            i = min(row)
            if i not in pivots:
                scale = row.pop(i)
                pivots[i] = {j: -v/scale for j, v in row.items()}
                break
            scale = row.pop(i)
            add_scaled(row, pivots[i], scale)
    free = [i for i in range(size) if i not in pivots]
    images = {i: {j: Fraction(1)} for j, i in enumerate(free)}
    for i in sorted(pivots, reverse=True):
        images[i] = {}
        for j, value in pivots[i].items():
            add_scaled(images[i], images[j], value)
    return len(free), [images[i] for i in range(size)]


class PathRow:
    """Path quotients from one start; multiplication appends the last arrow."""
    def __init__(self, windows, start, length):
        self.dims = [1]
        self.arrows = []
        for d in range(length):
            relations = []
            if d:
                for b, relation in product(range(self.dims[d-1]), windows[(start+d-1) % 4]):
                    vector = {}
                    for (x, y), value in relation.items():
                        for i, c in self.arrows[-1][4*b+x].items():
                            add_scaled(vector, {4*i+y: c}, value)
                    relations.append(vector)
            size, images = quotient_images(4*self.dims[d], relations)
            self.dims.append(size)
            self.arrows.append(images)

    def dimension(self, degree):
        return self.dims[degree] if degree >= 0 else 0


class CutExt:
    def __init__(self, first, second):
        self.windows = complete_windows(first, second)
        self.dual = orthogonal_windows([first, second])

    @lru_cache(maxsize=128)
    def frame(self, kind, start, length):
        return word_frame(self.windows if kind == 'path' else self.dual, start, length)

    @lru_cache(maxsize=256)
    def extend(self, kind, start, length, side, index, letter):
        """Multiply in B or its quadratic dual by a single letter."""
        _, words = self.frame(kind, start, length)
        word = words[index]
        if side == 'left':
            start -= 1
            word += letter * 4**length
        else:
            word = 4*word + letter
        larger, _ = self.frame(kind, start, length+1)
        return {j: v[word] for j, v in enumerate(larger) if word in v}

    def complex(self, row, offset, endpoint):
        """M at vertex a is row[offset+a]. Return cochain dimensions and sparse columns."""
        if not 0 <= endpoint < 4:
            raise ValueError("expected a cut endpoint in 0..3")
        bases = []
        for p in range(4):
            basis = []
            for a in range(min(endpoint, 3-p)+1):
                path, _ = self.frame('path', a, endpoint-a)
                dual, _ = self.frame('dual', a, p)
                basis.extend((a, j, k, b) for j, k, b in product(
                    range(len(path)), range(len(dual)), range(row.dimension(offset+a+p))))
            bases.append(basis)
        maps = []
        for p in range(3):
            indices = {word: i for i, word in enumerate(bases[p+1])}
            columns = []
            for a, j, k, b in bases[p]:
                column = {}
                for x in range(4):
                    if a:
                        head = self.extend('path', a, endpoint-a, 'left', j, x)
                        middle = self.extend('dual', a, p, 'left', k, x)
                        for t, c in head.items():
                            for u, v in middle.items():
                                add_scaled(column, {indices[a-1,t,u,b]: c*v}, 1)
                    if a+p < 3:
                        middle = self.extend('dual', a, p, 'right', k, x)
                        tail = row.arrows[offset+a+p][4*b+x]
                        for u, v in middle.items():
                            for t, c in tail.items():
                                add_scaled(column, {indices[a,j,u,t]: c*v}, (-1)**(p+1))
                columns.append(column)
            maps.append(columns)
        return list(map(len, bases)), maps
