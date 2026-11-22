"""Function-field arithmetic of the genus-three octic in characteristic 109."""

from flint import nmod_poly
from functools import cache
from rank_fifteen_octic import octic_frame
from rank_fourteen_function_field import evaluate_relations
from verify_fixed_geometric_quartic_locus import a, e, f
from verify_uniform_fixed_quartic_residuals import b
from verify_uniform_fixed_quartic_support import z, l

P = 109


def poly(v):
    return nmod_poly(list(v), P)


ZERO = poly([])
ONE = poly([1])
X = poly([0, 1])
F = X**8 + 14 * X**4 + 1


class Function:
    """(a+bY)/d over F_109(t), with a monic common denominator."""

    __slots__ = ("a", "b", "d", "inv")

    def __init__(self, a, b=ZERO, d=ONE):
        if not d:
            raise ZeroDivisionError()
        if not a and not b:
            self.a = self.b = ZERO
            self.d = ONE
        else:
            if d != ONE:
                g = d.gcd(a).gcd(b)
                if g.degree() > 0:
                    a, b, d = a // g, b // g, d // g
                lc = int(d.leading_coefficient())
                if lc != 1:
                    u = pow(lc, -1, P)
                    a, b, d = a * u, b * u, d * u
            self.a, self.b, self.d = a, b, d
        self.inv = None

    @classmethod
    def lift(cls, v):
        return v if isinstance(v, cls) else cls(poly([v]))

    def __bool__(self):
        return bool(self.a) or bool(self.b)

    def __neg__(self):
        return Function(-self.a, -self.b, self.d)

    def __add__(self, v):
        v = self.lift(v)
        if self.d == v.d:
            return Function(self.a + v.a, self.b + v.b, self.d)
        g = self.d.gcd(v.d)
        s = self.d // g
        t = v.d // g
        return Function(self.a * t + v.a * s, self.b * t + v.b * s, self.d * t)

    __radd__ = __add__

    def __sub__(self, v):
        return self + -self.lift(v)

    def __rsub__(self, v):
        return -self + v

    def __mul__(self, v):
        v = self.lift(v)
        return Function(
            self.a * v.a + self.b * v.b * F, self.a * v.b + self.b * v.a, self.d * v.d
        )

    __rmul__ = __mul__

    def norm_numerator(self):
        return self.a * self.a - self.b * self.b * F

    def inverse(self):
        if self.inv is None:
            self.inv = Function(
                self.a * self.d, -self.b * self.d, self.norm_numerator()
            )
        return self.inv

    def __truediv__(self, v):
        return self * self.lift(v).inverse()

    def __rtruediv__(self, v):
        return self.inverse() * v

    def __pow__(self, n):
        if n < 0:
            return self.inverse() ** (-n)
        if not n:
            return Function(ONE)
        v = self ** (n // 2)
        sq = v * v
        return sq * self if n % 2 else sq


@cache
def windows(color, j):
    fr = octic_frame(Function(X), Function(ZERO, ONE), Function.lift(j), color)
    return evaluate_relations(dict(zip((z, l, a, b, e, f), fr)), Function(ONE))
