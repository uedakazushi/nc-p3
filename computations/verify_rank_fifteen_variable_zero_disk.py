"""The first possible path obstruction on the whole small-trace, small-y disk at 3."""

from functools import cache
from itertools import product

from flint import fmpz_mod_poly_ctx, fmpz_poly, fq_default_ctx

from verify_klein_center_lift import central_equations, field_paths
from verify_rank_fifteen_cubic_zero_disk import UNCUT, pivot_rows, tensor_relations
from verify_rank_fifteen_parameter_orbits import pattern
from verify_rank_fifteen_quadratic_dual import sixth_insertion
from verify_rank_fourteen_all_ordinary import DIMENSIONS, saturated_projection

F9 = fq_default_ctx(modulus=fmpz_mod_poly_ctx(3)([1, 0, 1]), var="i")
BASE_POINT = 1 + F9.gen()
ONE = fmpz_poly([1])
FIXED_CENTER_WORDS = (("1313", "2222", "2323", "2332", "3131", "3223", "3232", "3322", "3333"),
                      ("1312", "2223", "2232", "2322", "2333", "3130", "3222", "3233", "3323", "3332"))


@cache
def residue(poly):
    result = F9.zero()
    for c in reversed(poly):
        result = result*BASE_POINT + c
    return result


class Rational:
    """A reduced fraction over Z[z], with evaluation at the chosen geometric residue."""

    __slots__ = ("n", "d", "value")

    def __init__(self, numerator, denominator=ONE):
        if not denominator:
            raise ZeroDivisionError
        if denominator != ONE:
            common = numerator.gcd(denominator)
            if denominator.leading_coefficient() < 0:
                common = -common
            numerator, denominator = numerator // common, denominator // common
        self.n, self.d, self.value = numerator, denominator, None

    @staticmethod
    @cache
    def integer(n):
        return Rational(fmpz_poly([n]))

    def __bool__(self):
        return bool(self.n)

    def __neg__(self):
        return Rational(-self.n, self.d)

    def __add__(self, other):
        if not isinstance(other, Rational):
            other = self.integer(other)
        if not self:
            return other
        if not other:
            return self
        if self.d == other.d:
            return Rational(self.n+other.n, self.d)
        return Rational(self.n*other.d + other.n*self.d, self.d*other.d)

    __radd__ = __add__

    def __mul__(self, other):
        if not isinstance(other, Rational):
            other = self.integer(other)
        if not self or not other:
            return self.integer(0)
        return Rational(self.n*other.n, self.d*other.d)

    __rmul__ = __mul__

    def inverse(self):
        return Rational(self.d, self.n)

    def specialize(self):
        if self.value is None:
            denominator = residue(tuple(map(int, self.d)))
            assert not denominator.is_zero()
            self.value = residue(tuple(map(int, self.n))) / denominator
        return self.value


class Jet:
    """Truncated H-series over Q(z); the least reliable input sets the precision."""

    __slots__ = ("cs",)
    precision = 7

    def __init__(self, coefficients):
        self.cs = tuple(coefficients)

    @classmethod
    def scalar(cls, coefficient):
        return cls([coefficient] + [Rational.integer(0)]*(cls.precision-1))

    @classmethod
    def constant(cls, n):
        return cls.scalar(Rational.integer(n))

    def __bool__(self):
        return any(self.cs)

    def __neg__(self):
        return Jet([-c for c in self.cs])

    def __add__(self, other):
        if not isinstance(other, Jet):
            other = self.constant(other)
        n = min(len(self.cs), len(other.cs), self.precision)
        return Jet([self.cs[j]+other.cs[j] for j in range(n)])

    __radd__ = __add__

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        if not isinstance(other, Jet):
            return Jet([c*other for c in self.cs])
        n = min(len(self.cs), len(other.cs), self.precision)
        return Jet([sum(self.cs[j]*other.cs[k-j] for j in range(k+1)) for k in range(n)])

    __rmul__ = __mul__

    def unit(self):
        return not self.cs[0].specialize().is_zero()

    def inverse(self):
        leading_inverse = self.cs[0].inverse()
        result = [leading_inverse]
        for k in range(1, min(self.precision, len(self.cs))):
            result.append(-leading_inverse * sum(self.cs[j]*result[k-j] for j in range(1, k+1)))
        return Jet(result)

    def divide_by_parameter(self):
        assert not self.cs[0]  # An identity in Q(z).
        return Jet([*self.cs[1:], Rational.integer(0)])


def relations():
    one = Jet.constant(1)
    half = Jet.scalar(Rational(ONE, fmpz_poly([2])))
    z = Jet.scalar(Rational(fmpz_poly([0, 1])))
    h = Jet([Rational.integer(j == 1) for j in range(Jet.precision)])
    return trace_relations(one, half, z, h, Jet.precision)


def trace_relations(one, half, z, h, precision):
    zi = z.inverse()
    p0, q0 = z+zi, z-zi
    K, n0 = (z*z+zi*zi)*half, (z*z-zi*zi)*half

    def root(square, initial):
        result = initial
        for _ in range(precision.bit_length()):
            result = (result + square*result.inverse())*half
        assert not (result*result-square)
        return result

    p, u = (root(p0*p0 + sign*2*h, p0) for sign in (1, -1))
    q, v = (root(q0*q0 + sign*2*h, q0) for sign in (1, -1))
    n2 = n0*n0
    n = root((n2 + root(n2*n2-4*h*h, n2))*half, n0)
    gamma = pattern(one, (K-1)*half, n*half, h*n.inverse()*half,
                    (u-p)*half*half, (u+p)*half*half, (v-q)*half*half, (v+q)*half*half)
    return tensor_relations(gamma, one*0)


def center_data(precision, jet=Jet, build_relations=relations, *, fixed_sections=True):
    jet.precision = precision
    windows = build_relations()
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

    basis, image = field_paths(windows, jet.constant(1), reduce)
    assert jet.precision == precision - 4
    centers = [[] for _ in range(4)]
    for h, dimension in enumerate((3, 2, 0, 0)):
        rows, width = central_equations(basis, image, h)
        pivots, loss = pivot_rows(rows, width, dimension, jet)
        assert not loss
        free = [j for j in range(width) if j not in pivots]
        if fixed_sections and h in (0, 1):
            # These free coordinates fix the same exact center sections when
            # the jet precision changes. All are in the last start's word basis.
            labels = tuple("".join(map(str, w)) for w in basis[3, 4, h])
            assert labels == FIXED_CENTER_WORDS[h]
            assert free == list(range(width-dimension, width))
        for col in (free[1:] if h == 0 else free):
            vector = {col: jet.constant(1)}
            vector.update({j: -row[col] for j, row in pivots.items() if col in row})
            offset = 0
            for p in range(4):
                words = basis[p, 4, h]
                centers[p].append((h, {w: vector[offset+j] for j, w in enumerate(words)
                                       if offset+j in vector}))
                offset += len(words)

    return windows, basis, image, centers


def verify():
    windows, basis, image, centers = center_data(7)
    for p, h in product(range(4), repeat=2):
        words = basis[p, 4, h]
        rows = [{j: center[w] for j, w in enumerate(words) if w in center}
                for g, center in centers[p] if g == h]
        _, _, loss = saturated_projection(rows, len(words), DIMENSIONS[4][h], Jet)
        assert loss == 0
        rows, width = sixth_insertion(windows, basis, image, p, h)
        assert width == 56
        _, _, loss = saturated_projection(rows, width, (20, 20, 22, 22)[h], Jet)
        assert loss == 2
    print("Variable trace at three: all sixth insertions, center ranks, and degree-four cuts verified.")


if __name__ == "__main__":
    verify()
