"""All four central cuts on the slope 59319/10000 and the indexed pencil geometry."""

from concurrent.futures import ProcessPoolExecutor

import sympy as sp
from flint import ctx, fmpq_poly

from verify_rank_fifteen_cubic_zero_disk import PRECISION, check_cut, tensor_relations
from verify_rank_fifteen_parameter_orbits import pattern
from verify_rank_fourteen_all_ordinary import gaussian_jets


def pencil_geometry():
    D, x, s = sp.symbols("D x s")
    K = 3*sp.I/4
    slope = sp.Rational(59319, 10000)
    assert sp.cancel((K*K+3)**3/(K*K-1)**2) == slope
    assert slope not in (0, 27) and slope**2-76*slope+1372 != 0
    branch = (s*(s*s-1)/(9*s*s-1))**2
    assert sp.cancel(((9*s*s+3)**3/(9*s*s-1)**2-27)/3**6-branch) == 0
    even_cubic = (1-x)**2-4*x**3/D
    shift_cubic = x*(x-1)**2-4/D
    assert sp.cancel(even_cubic.subs(x, s**-2)*s**6 - shift_cubic.subs(x, s**2)) == 0
    for cubic, j in ((even_cubic, D*(D-24)**3/(D-27)), (shift_cubic, 16*D**2/(D-27))):
        a, b, c, _ = sp.Poly(cubic, x).all_coeffs()
        disc = sp.discriminant(cubic, x)
        assert sp.cancel(disc - 16*(D-27)/D**2) == 0
        assert sp.cancel(256*(b*b-3*a*c)**3/(a*a*disc)-j) == 0


def verify_start(start):
    ctx.cap = PRECISION+5
    Jet = gaussian_jets(-1, PRECISION, residue_prime=3)
    one = Jet.constant(1)
    half = Jet(fmpq_poly([1])/2)
    h = Jet(fmpq_poly([0, 1]))
    z = Jet(fmpq_poly([1]), fmpq_poly([1]))
    zi = z.inverse()
    p0, q0 = z+zi, z-zi
    K, n0 = (z*z+zi*zi)*half, (z*z-zi*zi)*half

    def root(square, initial):
        value = initial
        for _ in range(PRECISION.bit_length()):
            value = (value+square*value.inverse())*half
        assert not (value*value-square)
        return value

    p, u = (root(p0*p0+sign*2*h, p0) for sign in (1, -1))
    q, v = (root(q0*q0+sign*2*h, q0) for sign in (1, -1))
    n2 = n0*n0
    n = root((n2+root(n2*n2-4*h*h, n2))*half, n0)
    gamma = pattern(one, (K-1)*half, n*half, h*n.inverse()*half,
                    (u-p)*half*half, (u+p)*half*half, (v-q)*half*half, (v+q)*half*half)
    losses, remaining = check_cut(start, tensor_relations(gamma, Jet.constant(0)), Jet)
    assert losses == [0, 0, 0, 2, 4, 2, 2, 2, 2, 2, 2, 2, 0] and remaining == 1
    return start


if __name__ == "__main__":
    pencil_geometry()
    with ProcessPoolExecutor(max_workers=2) as pool:
        for start in pool.map(verify_start, range(4)):
            print("Verified Gaussian trace through the thirteenth cut at start", start, flush=True)
