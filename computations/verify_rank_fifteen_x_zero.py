"""Certify x=0 at the ramified trace K=i*sqrt(3) on the small-y disk at three."""

from concurrent.futures import ProcessPoolExecutor
from functools import cache

from flint import ctx, fmpq_poly, fmpq_series

from verify_rank_fifteen_cubic_zero_disk import PRECISION, ZERO, check_cut, tensor_relations
from verify_rank_fifteen_parameter_orbits import pattern


class Jet:
    """Coefficients of 1,i,sqrt(3),i*sqrt(3) as truncated rational polynomials in H."""

    __slots__ = ("coefficients",)
    precision = PRECISION

    def __init__(self, coefficients):
        self.coefficients = coefficients

    @staticmethod
    @cache
    def constant(value):
        return Jet((fmpq_poly([value]), ZERO, ZERO, ZERO))

    def __bool__(self):
        return any(self.coefficients)

    def __add__(self, other):
        if not isinstance(other, Jet):
            return Jet((self.coefficients[0] + other, *self.coefficients[1:]))
        return Jet(tuple(a+b for a, b in zip(self.coefficients, other.coefficients)))

    __radd__ = __add__

    def __neg__(self):
        return Jet(tuple(-a for a in self.coefficients))

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    @classmethod
    def gaussian_product(cls, a, b, c, d):
        n = cls.precision
        return a.mul_low(c, n) - b.mul_low(d, n), a.mul_low(d, n) + b.mul_low(c, n)

    def __mul__(self, other):
        if not isinstance(other, Jet):
            return Jet(tuple(c*other for c in self.coefficients))
        a, b, c, d = self.coefficients
        e, f, g, h = other.coefficients
        real, imag = self.gaussian_product(a, b, e, f)
        real_d, imag_d = self.gaussian_product(c, d, g, h)
        mixed, mixed_i = self.gaussian_product(a+c, b+d, e+g, f+h)
        return Jet((real+3*real_d, imag+3*imag_d,
                    mixed-real-real_d, mixed_i-imag-imag_d))

    __rmul__ = __mul__

    def unit(self):
        constants = [c[0] for c in self.coefficients]
        assert all(int(x.denominator) % 3 for x in constants)
        # sqrt(3) vanishes in the residue field, where 1,i are independent.
        return any(int(x.numerator) % 3 for x in constants[:2])

    def inverse(self):
        assert self.unit()
        a, b, c, d = self.coefficients
        real, imag = self.gaussian_product(a, b, a, b)
        real_d, imag_d = self.gaussian_product(c, d, c, d)
        real, imag = real-3*real_d, imag-3*imag_d
        n = self.precision
        norm = real.mul_low(real, n) + imag.mul_low(imag, n)
        inverse = fmpq_poly(list(1 / fmpq_series(norm, prec=n)))
        return Jet(tuple(v.mul_low(inverse, n) for v in
                         (*self.gaussian_product(a, b, real, -imag),
                          *self.gaussian_product(-c, -d, real, -imag))))

    def divide_by_parameter(self):
        assert all(not c[0] for c in self.coefficients)
        return Jet(tuple(c.right_shift(1) for c in self.coefficients))


def check_start(start):
    Jet.precision = PRECISION
    ctx.cap = PRECISION + 5
    one = Jet.constant(1)
    unit = fmpq_poly([1])
    i = Jet((ZERO, unit, ZERO, ZERO))
    d = Jet((ZERO, ZERO, unit, ZERO))
    h = Jet((fmpq_poly([0, 1]), ZERO, ZERO, ZERO))
    half = Jet((unit/2, ZERO, ZERO, ZERO))
    K = i*d

    def root(value, initial):
        ans = initial
        for _ in range(PRECISION.bit_length()):
            ans = (ans + value*ans.inverse())*half
        assert not (ans*ans - value)
        return ans

    p, u = (root(2*(K+1+sign*h), d+i) for sign in (1, -1))
    q, v = (root(2*(K-1+sign*h), -(1+K)) for sign in (1, -1))
    n = root((-4-root(16-4*h*h, 4*one))*half, -2*i)
    gamma = pattern(one, (K-1)*half, n*half, h*n.inverse()*half,
                    (u-p)*half*half, (u+p)*half*half, (v-q)*half*half, (v+q)*half*half)
    return check_cut(start, tensor_relations(gamma, Jet.constant(0)), Jet)


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=3) as pool:
        for p, (losses, precision) in enumerate(pool.map(check_start, range(4))):
            print("Start", p, "losses", losses, "remaining precision", precision, flush=True)
