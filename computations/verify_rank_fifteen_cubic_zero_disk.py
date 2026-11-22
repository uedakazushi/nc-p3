"""Saturate the four central cuts on the K=3, H=0 disk above the prime three."""

from concurrent.futures import ProcessPoolExecutor
from functools import cache
from itertools import product

from flint import ctx, fmpq_poly, fmpq_series

from verify_klein_center_lift import central_equations, field_paths
from verify_rank_fifteen_parameter_orbits import dot, pattern
from verify_rank_fourteen_all_ordinary import DIMENSIONS, saturated_projection
from verify_rank_fourteen_central_cut import quotient_layers

PRECISION = 25
ZERO = fmpq_poly([])
UNCUT = ((1, 0, 0, 0), (1,)*4, (2, 2, 3, 3), (5,)*4, (9, 10, 8, 8), (14,)*4)


class Jet:
    """Truncated a(H)+sqrt(2)b(H), integral at the inert prime over three."""

    __slots__ = ("a", "b")
    precision = PRECISION

    def __init__(self, a, b=ZERO):
        self.a, self.b = a, b

    @staticmethod
    @cache
    def constant(value):
        return Jet(fmpq_poly([value]))

    def __bool__(self):
        return bool(self.a or self.b)

    def __neg__(self):
        return Jet(-self.a, -self.b)

    def __add__(self, other):
        if not isinstance(other, Jet):
            other = self.constant(other)
        return Jet(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        if not isinstance(other, Jet):
            return Jet(self.a*other, self.b*other)
        n = self.precision
        return Jet(
            self.a.mul_low(other.a, n) + 2*self.b.mul_low(other.b, n),
            self.a.mul_low(other.b, n) + self.b.mul_low(other.a, n),
        )

    __rmul__ = __mul__

    def unit(self):
        a, b = self.a[0], self.b[0]
        assert int(a.denominator) % 3 and int(b.denominator) % 3
        return bool(int(a.numerator) % 3 or int(b.numerator) % 3)

    def inverse(self):
        assert self.unit()
        n = self.precision
        norm = self.a.mul_low(self.a, n) - 2*self.b.mul_low(self.b, n)
        reciprocal = fmpq_poly(list(1 / fmpq_series(norm, prec=n)))
        return Jet(reciprocal.mul_low(self.a, n), -reciprocal.mul_low(self.b, n))

    def divide_by_parameter(self):
        # Check exact vanishing in Q(sqrt(2)), including all multiples of three.
        assert not self.a[0] and not self.b[0]
        return Jet(self.a.right_shift(1), self.b.right_shift(1))


def windows():
    one = Jet.constant(1)
    half = Jet(fmpq_poly([1])/2)
    h = Jet(fmpq_poly([0, 1]))
    sqrt_two = Jet(ZERO, fmpq_poly([1]))

    def root(value, initial):
        approximation = initial
        for _ in range(PRECISION.bit_length()):
            approximation = (approximation + value*approximation.inverse())*half
        assert not (approximation*approximation - value)
        return approximation

    p, u = (root(8 + sign*2*h, 2*sqrt_two) for sign in (1, -1))
    q, v = (root(4 + sign*2*h, 2*one) for sign in (1, -1))
    n = root((8 + root(64 + (-4)*h*h, 8*one))*half, 2*sqrt_two)
    gamma = pattern(one, one, n*half, h*n.inverse()*half, (u-p)*half*half,
                    (u+p)*half*half, (v-q)*half*half, (v+q)*half*half)

    return tensor_relations(gamma, Jet.constant(0))


def tensor_relations(gamma, zero):
    def component(a, b, c, d):
        if a ^ c != b ^ d:
            return zero
        return sum(
            gamma[(a ^ c) + 4*j] * (-1)**(dot(j, c) + dot((0, 3, 2, 1)[j], d))
            for j in range(4)
        )

    result = []
    for start in range(4):
        side = []
        for color, k in product(range(4), repeat=2):
            row = {}
            for j in range(4):
                word = (j, j ^ color, k, k ^ color)
                value = component(*(word[-start:] + word[:-start]))
                if value:
                    row[j, j ^ color] = value
            if row:
                side.append((color, row))
        result.append(side)
    return result


def pivot_rows(rows, width, target, jet):
    free, projection, loss = saturated_projection(rows, width, target, jet)
    pivots = {}
    for j in range(width):
        if j not in free:
            pivots[j] = {free[k]: -v for k, v in projection[j].items()}
            pivots[j][j] = jet.constant(1)
    return pivots, loss


def four_sections(relations, jet):
    # field_paths visits all sixteen (start,color) blocks before each next length.
    order = iter(product(range(1, 6), range(4), range(4)))
    losses = []

    def reduce(rows):
        n, p, h = next(order)
        pivots, loss = pivot_rows(rows, sum(UNCUT[n-1]), UNCUT[n][h], jet)
        losses.append(loss)
        if p == h == 3:
            assert max(losses) == (0, 0, 0, 2, 2)[n-1]
            jet.precision -= max(losses)
            losses.clear()
        return pivots

    one = jet.constant(1)
    basis, image = field_paths(relations, one, reduce)
    assert next(order, None) is None
    sections = [[] for _ in range(4)]
    for h, dimension in enumerate((3, 2, 0, 0)):
        rows, width = central_equations(basis, image, h)
        pivots, loss = pivot_rows(rows, width, dimension, jet)
        assert loss == 0
        free = [j for j in range(width) if j not in pivots]
        for col in (free[1:] if h == 0 else free):
            vector = {col: one}
            vector.update({j: -row[col] for j, row in pivots.items() if col in row})
            offset = 0
            for p in range(4):
                words = basis[p, 4, h]
                sections[p].append((h, {w: vector[offset+k] for k, w in enumerate(words)
                                        if offset+k in vector}))
                offset += len(words)
    return sections


def check_start(start):
    Jet.precision = PRECISION
    ctx.cap = PRECISION + 5
    return check_cut(start, windows(), Jet)


def check_cut(start, relations, jet):
    sections = four_sections(relations, jet)
    losses = []

    def project(n, h, rows, size):
        free, projection, loss = saturated_projection(rows, size, DIMENSIONS[n][h], jet)
        losses.append(loss)
        return free, projection

    losses_by_length = []
    for _ in quotient_layers(relations, sections, start, 13, jet.constant(1), project):
        losses_by_length.append(max(losses))
        jet.precision -= max(losses)
        assert jet.precision > 0
        losses.clear()
    return losses_by_length, jet.precision


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=2) as pool:
        for p, (losses, precision) in enumerate(pool.map(check_start, range(4))):
            print("Start", p, "losses", losses, "remaining precision", precision, flush=True)
