"""Exact F3(c)[u] arithmetic on the residual sextic, with pole tracking."""

from flint import nmod_poly

ZERO = nmod_poly([], 3)
ONE = nmod_poly([1], 3)
C = nmod_poly([0, 1], 3)

# Write t=u^2. The entries are the coefficients of 1,t,t^2 in t^3.
REDUCTION = (-(C**3) + C**2 + C - 1, -(C**4) + C**3 - C, C**2 + 1)
ZERO_TRIPLE = (ZERO, ZERO, ZERO)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b, strict=True))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b, strict=True))


def scale(a, coefficient):
    return tuple(x * coefficient for x in a)


def times_t(a):
    x, y, z = a
    return z * REDUCTION[0], x + z * REDUCTION[1], y + z * REDUCTION[2]


def multiply(a, b):
    high = a[2] * b[2]
    middle = a[1] * b[2] + a[2] * b[1] + high * REDUCTION[2]
    return (
        a[0] * b[0] + middle * REDUCTION[0],
        a[0] * b[1] + a[1] * b[0] + middle * REDUCTION[1] + high * REDUCTION[0],
        a[0] * b[2]
        + a[1] * b[1]
        + a[2] * b[0]
        + middle * REDUCTION[2]
        + high * REDUCTION[1],
    )


def adjugate_column(a):
    """First column of the adjugate of multiplication by a, and its determinant."""
    second = times_t(a)
    third = times_t(second)
    column = (
        second[1] * third[2] - second[2] * third[1],
        a[2] * third[1] - a[1] * third[2],
        a[1] * second[2] - a[2] * second[1],
    )
    det = a[0] * column[0] + second[0] * column[1] + third[0] * column[2]
    return column, det


class Curve:
    """(even(t)+u*odd(t))/den, with a common monic polynomial denominator."""

    __slots__ = ("even", "odd", "den")
    divisor = ONE

    def __init__(self, even, odd=ZERO_TRIPLE, den=ONE):
        if not den:
            raise ZeroDivisionError
        if den != ONE:
            if not any(even + odd):
                even, odd, den = ZERO_TRIPLE, ZERO_TRIPLE, ONE
            else:
                common = den
                for coeff in even + odd:
                    if common.degree() <= 0:
                        break
                    if coeff:
                        common = common.gcd(coeff)
                if common.degree() > 0:
                    even = tuple(z // common for z in even)
                    odd = tuple(z // common for z in odd)
                    den = den // common
                if int(den.leading_coefficient()) == 2:
                    even, odd, den = scale(even, -1), scale(odd, -1), -den
        self.even, self.odd, self.den = even, odd, den

    @classmethod
    def constant(cls, x):
        if isinstance(x, cls):
            return x
        poly = x if isinstance(x, nmod_poly) else nmod_poly([int(x)], 3)
        return cls((poly, ZERO, ZERO))

    def __bool__(self):
        return any(self.even + self.odd)

    def __neg__(self):
        return Curve(scale(self.even, -1), scale(self.odd, -1), self.den)

    def __add__(self, x):
        x = self.constant(x)
        if not self:
            return x
        if not x:
            return self
        if self.den == x.den:
            return Curve(add(self.even, x.even), add(self.odd, x.odd), self.den)
        gcd = self.den.gcd(x.den)
        left = x.den // gcd
        right = self.den // gcd
        return Curve(
            add(scale(self.even, left), scale(x.even, right)),
            add(scale(self.odd, left), scale(x.odd, right)),
            self.den * left,
        )

    __radd__ = __add__

    def __sub__(self, x):
        return self + (-self.constant(x))

    def __rsub__(self, x):
        return (-self) + x

    def __mul__(self, x):
        x = self.constant(x)
        if not self or not x:
            return Curve(ZERO_TRIPLE)
        even = multiply(self.even, x.even)
        odd = multiply(self.odd, x.odd)
        cross = multiply(add(self.even, self.odd), add(x.even, x.odd))
        return Curve(
            add(even, times_t(odd)),
            subtract(subtract(cross, even), odd),
            self.den * x.den,
        )

    __rmul__ = __mul__

    def inverse(self):
        norm = subtract(
            multiply(self.even, self.even), times_t(multiply(self.odd, self.odd))
        )
        adj, den = adjugate_column(norm)
        if not den:
            raise ZeroDivisionError
        p = den
        while p.degree() > 0:
            gcd = p.gcd(Curve.divisor)
            if gcd.degree() <= 0:
                for factor, _ in p.factor_squarefree()[1]:
                    Curve.divisor *= factor * int(factor.leading_coefficient())
                break
            p = p // gcd
        return Curve(
            scale(multiply(self.even, adj), self.den),
            scale(multiply(self.odd, adj), -self.den),
            den,
        )

    def __truediv__(self, x):
        return self * self.constant(x).inverse()

    def __rtruediv__(self, x):
        return self.inverse() * x

    def __pow__(self, k):
        if k < 0:
            return self.inverse() ** (-k)
        if k == 0:
            return self.constant(1)
        if k == 1:
            return self
        half = self ** (k // 2)
        square = half * half
        return square * self if k % 2 else square

    def pivot_degree(self):
        return max(p.degree() for p in self.even + self.odd) + self.den.degree()
