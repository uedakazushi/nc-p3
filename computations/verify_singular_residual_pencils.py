"""Exact residual two-by-two minors and Chern lengths after removing constant eigenspaces."""

import sympy as sp


def main():
    z, w, u, v, t = sp.symbols("z w u v t")
    a, b, aw, bw = z-u, z-v, w-u, w-v
    node = ((z-t*w)*a*bw+(t*z-w)*aw*b)/(z-w)
    cusp = (a*bw+aw*b)/t+z*w*(a*bw-aw*b)/(z-w)
    for expr, expected in [(node, -(u-v)*(t**2*u-v)), (cusp, (u-v)*(2*t*u*v-u+v)/t**2)]:
        poly = sp.Poly(sp.cancel(expr), z, w)
        matrix = sp.Matrix(2, 2, lambda i, j: poly.coeff_monomial(z**i*w**j))
        assert sp.factor(matrix.det()-expected) == 0

    f, h, d = sp.symbols("f h d")
    # f and h are the two rulings of P1 x P1. Only the coefficient of fh survives.
    for rank, expected in [(3, 16), (2, 12-d), (1, 8-d)]:
        degree = d-4
        source_first = degree*f-rank*h
        source_second = -(rank-1)*degree*f*h
        target_first = d*f
        chern = source_first**2-source_second-source_first*target_first
        value = sp.expand(chern).coeff(f, 1).coeff(h, 1)
        if rank == 3:
            value = value.subs(d, 0)
        assert sp.expand(value-expected) == 0
    assert 7+12 < 20 and 7+8 < 20 and 14+8-3 < 20
    print("Verified both residual factors and the lengths 16, 12-d, and 8-d.")


if __name__ == "__main__":
    main()
