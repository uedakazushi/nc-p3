"""Saturated four-cut bounds at all places of trace ±1 above 109."""

from flint import fmpq_poly, fmpq_series, nmod_poly, ctx
from functools import cache
from itertools import product
from verify_rank_fourteen_all_ordinary import saturated_projection, DIMENSIONS
from verify_rank_fourteen_central_cut import quotient_layers
from rank_fifteen_octic_integral import integral_presentation
from concurrent.futures import ProcessPoolExecutor

P = 109
ZERO = fmpq_poly([])
ONE = fmpq_poly([1])
PARAMETER = fmpq_poly([0, 1])
RESIDUE_PARAMETER = nmod_poly([0, 1], P)


@cache
def residue(c):
    n, d = int(c.numerator), int(c.denominator)
    assert d % P
    return n * pow(d, -1, P) % P


def jets(N):
    """Truncated Q(i,√3)-series and one fast residue for rejecting nonunits."""
    ctx.cap = N + 5

    def multiply_low(a, b):
        return a.mul_low(b, Jet.precision) if a and b else ZERO

    def gaussian_product(a, b, c, d):
        return multiply_low(a, c) - multiply_low(b, d), multiply_low(
            a, d
        ) + multiply_low(b, c)

    class Jet:
        __slots__ = ("coefficients", "residue_series")
        precision = N

        def __init__(self, cs, mod):
            self.coefficients = cs
            self.residue_series = mod

        @classmethod
        def constant(cls, v):
            return Jet((fmpq_poly([v]), ZERO, ZERO, ZERO), nmod_poly([v], P))

        def __bool__(self):
            return any(self.coefficients)

        def __neg__(self):
            return Jet(tuple(-c for c in self.coefficients), -self.residue_series)

        def __add__(self, w):
            if not isinstance(w, Jet):
                if not w:
                    return self
                return Jet(
                    (self.coefficients[0] + w, *self.coefficients[1:]),
                    self.residue_series + w,
                )
            return Jet(
                tuple(a + b for a, b in zip(self.coefficients, w.coefficients)),
                self.residue_series + w.residue_series,
            )

        __radd__ = __add__

        def __sub__(self, w):
            return self + -w

        def __rsub__(self, w):
            return -self + w

        def __mul__(self, w):
            if not isinstance(w, Jet):
                return Jet(
                    tuple(a * w for a in self.coefficients), self.residue_series * w
                )
            a, b, c, d = self.coefficients
            e, f, g, h = w.coefficients
            A, B = gaussian_product(a, b, e, f)
            C, D = gaussian_product(c, d, g, h)
            if not c and not d:
                E, F = gaussian_product(a, b, g, h)
            elif not g and not h:
                E, F = gaussian_product(c, d, e, f)
            else:
                E, F = gaussian_product(a + c, b + d, e + g, f + h)
                E, F = E - A - C, F - B - D
            return Jet(
                (A + 3 * C, B + 3 * D, E, F),
                self.residue_series.mul_low(w.residue_series, Jet.precision),
            )

        __rmul__ = __mul__

        def inverse(self):
            assert self.unit()
            a, b, c, d = self.coefficients
            R, S = gaussian_product(a, b, a, b)
            T, U = gaussian_product(c, d, c, d)
            R, S = R - 3 * T, S - 3 * U
            norm = multiply_low(R, R) + multiply_low(S, S)
            inv = fmpq_poly(list(1 / fmpq_series(norm, prec=Jet.precision)))
            R, S = multiply_low(R, inv), -multiply_low(S, inv)
            p, q = gaussian_product(a, b, R, S)
            r, s = gaussian_product(-c, -d, R, S)
            return Jet(
                (p, q, r, s), self.residue_series.inverse_series_trunc(Jet.precision)
            )

        def unit(self):
            # Inversion must be integral at all four primes above 109.
            if not self.residue_series[0]:
                return False
            vals = [residue(poly[0]) for poly in self.coefficients]
            a, b, c, d = vals
            return all(
                (a + ri * b + rd * (c + ri * d)) % 109
                for ri in (33, -33)
                for rd in (49, -49)
            )

        def divide_by_parameter(self):
            assert (
                all(not a[0] for a in self.coefficients) and not self.residue_series[0]
            )
            return Jet(
                tuple(a.right_shift(1) for a in self.coefficients),
                self.residue_series.right_shift(1),
            )

    return (
        Jet,
        Jet((ZERO, ONE, ZERO, ZERO), nmod_poly([33], P)),
        Jet((ZERO, ZERO, ONE, ZERO), nmod_poly([49], P)),
        Jet((PARAMETER, ZERO, ZERO, ZERO), RESIDUE_PARAMETER),
    )


def expand(data, Jet, I, D, u, k, sgna, sgnb):
    """Expand the fixed integral equations at Y=0 and remove common powers of Y."""
    U, V = (sgna * D, sgnb * I) if k == 1 else (sgna * I, sgnb * D)
    alpha = (V - 1) * (U - 1).inverse()
    t = alpha

    def power(a, n):
        if not n:
            return Jet.constant(1)
        s = power(a, n // 2)
        return s * s * a if n % 2 else s * s

    for _ in range(Jet.precision.bit_length()):
        f = power(t, 8) + 14 * power(t, 4) + 1 - u * u
        if not f:
            break
        t = t - f * (8 * power(t, 7) + 56 * power(t, 3)).inverse()
    assert not (power(t, 8) + 14 * power(t, 4) + 1 - u * u)

    @cache
    def ev(poly):
        v = Jet.constant(0)
        for c in reversed(poly):
            v = v * t + c
        return v

    def entry(coords):
        a, b, c, d, e, f, g, h = [ev(tuple(q)) for q in coords]
        return a + I * b + D * (e + I * f) + u * (c + I * d + D * (g + I * h))

    out = [
        [
            [(h, {w: entry(v) for w, v in row.items()}) for h, row in side]
            for side in family
        ]
        for family in data
    ]

    def normalize(rows):
        # Remove a power only after exact vanishing, also for a four-position center.
        d = 0
        while not any(v.unit() for row in rows for v in row.values()):
            d += 1
            assert d < Jet.precision
            for row in rows:
                for w, v in row.items():
                    row[w] = v.divide_by_parameter()
        for row in rows:
            for w in list(row):
                if not row[w]:
                    del row[w]
        return d

    ls = [normalize([row]) for side in out[0] for h, row in side]
    for column in range(4):
        ls.append(normalize([side[column][1] for side in out[1]]))
    return out, max(ls)


def verify_boundary(job):
    color, k, sign = job
    N = 24
    data = integral_presentation(color)
    Jet, I, D, u = jets(N)
    (ws, cs), initial_loss = expand(data, Jet, I, D, u, k, 1, sign)
    for start in range(4):
        Jet.precision = N - initial_loss
        losses = []
        history = []

        def proj(n, h, rows, size):
            free, project, m = saturated_projection(rows, size, DIMENSIONS[n][h], Jet)
            losses.append(m)
            return free, project

        for n in quotient_layers(ws, cs, start, 13, Jet.constant(1), proj):
            history.append(max(losses))
            Jet.precision -= max(losses)
            assert Jet.precision > 0
            losses.clear()
        print(
            "COMPLETE",
            color,
            k,
            sign,
            start,
            initial_loss,
            history,
            Jet.precision,
            flush=True,
        )


def verify_color(color):
    for k, sign in product((1, -1), repeat=2):
        verify_boundary((color, k, sign))


if __name__ == "__main__":
    assert (33**2 + 1) % P == 0 and (49**2 - 3) % P == 0
    with ProcessPoolExecutor(max_workers=2) as pool:
        list(pool.map(verify_color, range(2)))
