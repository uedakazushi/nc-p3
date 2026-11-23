"""An integral 32-digit certificate on the outer pole branch H=m w^2."""

from flint import fmpz_poly, nmod_poly
import sympy as sp

from verify_rank_fifteen_integral_annulus import Recorded
from verify_rank_fifteen_pole_fibers import RecordingExact, cuts, tables
from verify_rank_fifteen_pole_slope_annulus import pattern, tensor_relations


def outer_relations(jet, half, w, u):
    one = jet.constant(1)
    h = u*w

    def root(square, initial):
        value = initial
        for _ in range(jet.precision.bit_length()):
            value = (value+square*value.inverse())*half
        assert not (value*value-square)
        return value

    denominator = (one-2*h*h).inverse()
    e = root(one-w*w, one)
    f = root(one-w*w-h*h, one)
    across = w*root(one-u*u, one)
    n, b = (2*h*denominator*(e+sign*w) for sign in (-1, 1))
    r2, t = 8*h*h*denominator*denominator, (one-2*w*w)*half
    assert not (n*b-r2*t) and not (n*n+b*b-r2)
    p, u_root = (root(4*denominator+sign*2*t*r2, jet.constant(2)) for sign in (1, -1))
    q, v = 4*h*denominator*f, 4*h*denominator*across
    gamma = pattern(one, 2*h*h*denominator, n*half, b*half,
                    (u_root-p)*half*half, (u_root+p)*half*half, (v-q)*half*half, (v+q)*half*half)
    return tensor_relations(gamma, one*0)


def relations(jet, half, w, m):
    return outer_relations(jet, half, w, m*w)


def verify():
    s, m = sp.symbols("s m")
    d, t = 1/(1-2*m*m*s*s), 1-2*s
    x, y = (1+2*m*m*s*s*d*d)*t**2, m**4*s**4*d**4*t**6
    leading = sp.cancel((x-1)/s)
    slope = sp.cancel(y/(x-1)**4)
    assert leading.subs(s, 0) == -4 and slope.subs(s, 0) == m**4/256
    assert sp.diff(slope, m).subs({s:0, m:1}) == sp.Rational(1,64)
    t = nmod_poly([0, 1], 3)
    baseline = sum((t**j for j in range(1, 9)), 0*t)
    other = t*t+t+2
    assert (t**9-t) == (t-1)*baseline and not (t**4+1)%other
    for modulus in ((-1, 1), (2, 1, 1)):
        class Outer(RecordingExact):
            pass
        Outer.modulus = fmpz_poly(list(modulus))
        Recorded.divisor = 1+0*t
        Recorded.record.cache_clear()
        products, _ = tables(
            Outer, 32, exact=True, build=relations, uncut_losses=(0, 3, 3, 4, 4),
            center_losses=(0, 2, 2, 2), fifth_loss=3,
        )
        cuts(Outer, products, losses=(3, 1, 2, 1, 3, 0, 2, 0))
        assert Outer.precision == 1
        if len(modulus) == 2:
            assert Recorded.divisor == baseline
        else:
            assert other.gcd(Recorded.divisor).degree() == 0
        print("Verified the exact outer chart at", modulus, "through degree thirteen.", flush=True)


if __name__ == "__main__":
    verify()
