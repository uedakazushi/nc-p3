"""The elliptic equality boundary at negative x and a 21-digit diagonal neighborhood."""

from collections import Counter

from flint import ctx, fmpq_poly, fq_default_ctx, fq_default_poly_ctx, nmod_poly

from verify_rank_fifteen_pole_fibers import RecordingPole, cuts
from verify_rank_fifteen_pole_slope_annulus import initial_tables, pattern, tensor_relations
from verify_rank_fifteen_unit_pole_depth import conic_profile
from verify_rank_fifteen_unit_residues import Scalar, graded_profile
from verify_rank_fifteen_vertical_all_levels import (
    EllipticUnit, ONE, T, ZERO, homogeneous_pullback, pivot_divisor,
    residue_polynomial, transform,
)
from verify_rank_fourteen_all_ordinary import gaussian_jets


def boundary_windows():
    field = EllipticUnit
    t = field.lift(T)
    i = field((ZERO, ONE, ZERO, ZERO))
    w = field((ZERO, ZERO, ONE, ZERO))
    a, b = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    g = 2*(1-i)*w/(1+t*t)
    p, u, q, v = 2*(a-b), 2*(a+b), g*(1-i), g*(1+i)
    gamma = pattern(field.constant(1), field.constant(0), -i*g/2, g/2,
                    (u-p)/4, (u+p)/4, (v-q)/4, (v+q)/4)
    return tensor_relations(gamma, field.constant(0))


def bad_slopes(divisor):
    num, den = 4*(T-T**3)**2, (1+T*T)**4
    assert num.gcd(den) == ONE
    core = divisor // divisor.gcd(num*den*(num-den))
    for a, b in ((-T, ONE), (ONE, T), (1+T, 1-T)):
        assert transform(num, a, b)*b**(8-num.degree())*den == transform(den, a, b)*num
    while core.degree() > 0:
        degree = core.degree()
        reflected = nmod_poly([int(c)*(-1)**(j%2) for j, c in enumerate(core)], 3)
        core = core.gcd(reflected)
        core = core.gcd(core.reverse())
        core = core.gcd(homogeneous_pullback(core, 1+T, 1-T))
        if core.degree() == degree:
            break

    fibers = {}
    for modulus, exponent in core.factor()[1]:
        assert exponent == 1
        gcd, inverse, _ = den.xgcd(modulus)
        assert gcd == ONE
        polynomial = residue_polynomial(num*inverse % modulus, modulus)
        fibers[polynomial] = fibers.get(polynomial, ONE)*modulus
    for polynomial, preimage in fibers.items():
        pullback = transform(nmod_poly(list(polynomial), 3), num, den)
        pullback *= pow(int(pullback.leading_coefficient()), -1, 3)
        assert pullback == preimage
    return core.degree(), tuple(sorted(fibers))


def negative_x_relations(field, half, w, g, *, imaginary_unit=None):
    one = field.constant(1)
    i = field(fmpq_poly([]), fmpq_poly([1])) if imaginary_unit is None else imaginary_unit
    s = i+w
    scale = (2*i+w).inverse()
    a, b = (1-s*s)*scale, 2*s*scale
    product = a*b

    def root(square, initial):
        value = initial
        for _ in range(field.precision.bit_length()):
            value = (value+square*value.inverse())*half
        assert not (value*value-square)
        return value

    k = root(one+w*w*g*g, one)
    p, u = (root(2*(k+1+sign*g*g*product), 2*one) for sign in (1, -1))
    q, v = (g*root(2*w*w*(k+1).inverse()+sign*2*product, one-sign*i)
            for sign in (1, -1))
    gamma = pattern(one, (k-1)*half, g*a*half, g*b*half,
                    (u-p)*half*half, (u+p)*half*half, (v-q)*half*half, (v+q)*half*half)
    return tensor_relations(gamma, 0*one)


def diagonal():
    precision = 21
    ctx.cap = precision+2
    field = gaussian_jets(-1, precision, residue_prime=3)
    orders = []
    field.record_pivot_orders = staticmethod(orders.append)
    t = field(fmpq_poly([0, 1]))
    products, frame = initial_tables(field, precision, field(fmpq_poly([1])/2), t, t,
                        build=negative_x_relations, uncut_losses=(0, 1, 1, 3, 2),
                        center_losses=(0, 1, 2, 2), fifth_loss=2)
    cuts(field, products, losses=(4, 1, 1, 1, 1, 0, 1, 0))
    assert field.precision == 1 and len(orders) == 308
    return frame, orders


def check_slope(polynomial):
    degree = len(polynomial)-1
    for multiplier in (1, 2, 4, 8, 16):
        extension = degree*multiplier
        if extension % 2:
            continue
        field = fq_default_ctx(3, extension, "a")
        ring = fq_default_poly_ctx(field)
        one = field.one()
        i = (-one).sqrt()
        slope = ring(list(polynomial)).roots()[0][0]
        for g, _ in (ring.gen()**4+16*slope).roots():
            h = -i*g*g
            if not all((4+sign*2*h).is_square() for sign in (1, -1)):
                continue
            p, u = (4+2*h).sqrt(), (4-2*h).sqrt()
            q, v = (1-i)*g, (1+i)*g
            assert g**4 == -16*slope and q*q == 2*h and v*v == -2*h
            gamma = pattern(one, field.zero(), -i*g/2, g/2,
                            (u-p)/4, (u+p)/4, (v-q)/4, (v+q)/4)
            entries = tuple(Scalar(field.zero() if c == 0 else c) for c in gamma)
            windows = tensor_relations(entries, Scalar(field.zero()))
            for pair in (((0, 1, 0), (0, 0, 1)), ((1, 0, 0), (0, 1, 0)),
                         ((1, 0, 1), (1, 1, 0))):
                try:
                    valid = graded_profile(windows, Scalar(one), pair)
                except AssertionError:
                    continue
                if valid:
                    print("Verified residual slope", polynomial, "in degree", extension, flush=True)
                    return
    raise AssertionError(f"No central quotient certificate for {polynomial}")


def verify():
    divisor = pivot_divisor(boundary_windows)
    assert divisor.degree() == 118054
    degree, exceptions = bad_slopes(divisor)
    assert degree == 840
    assert Counter(len(p)-1 for p in exceptions) == {1: 1, 2: 2, 3: 6, 4: 4, 5: 2, 6: 6, 8: 1, 12: 1}
    for polynomial in exceptions:
        check_slope(polynomial)
    frame, middle = diagonal()
    reference_frame, edge = conic_profile(RecordingPole)
    assert frame == reference_frame
    assert middle[:48] == edge[:48]
    assert sum(a != b for a, b in zip(middle[48:64], edge[48:64], strict=True)) == 8
    print("Verified the full negative-x equality boundary and the Gaussian diagonal tube.")


if __name__ == "__main__":
    verify()
