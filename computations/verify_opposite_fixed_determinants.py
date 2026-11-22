"""Signed sections with opposite block determinants on the positive fixed elliptic curve."""

import sympy as sp

from flint_rational_function import rational_function_field
from two_window_ext import PathRow
from verify_all_zero_pairing_defects import relation_bases
from verify_klein_elliptic_completion import coeffs, gram_twist


def section(r, l, sign):
    z = r * r
    a = r * (l * r * (r + 1) + sign * (r - 1))
    b = r + 1 + sign * l * r * (r - 1)
    e = r * (r + 1 - sign * l * r * (r - 1))
    f = l * r * (r + 1) - sign * (r - 1)
    return sp.Matrix([[a, b], [z * b, a]]), sp.Matrix([[e, f], [z * f, e]])


def verify_section(sign, r, l):
    z = r * r
    A, B = section(r, l, sign)
    a, b = A[0, :]
    e, f = B[0, :]
    x = a * a - z * f * f
    y = e * e - z * b * b
    assert sp.expand(x + y) == 0
    assert sp.expand(x - 4 * sign * l * r**3 * (z - 1)) == 0
    assert (
        sp.expand(A.det() - 4 * r**3 * (l * l * z - 1))
        == sp.expand(A.det() + B.det())
        == 0
    )
    co, (p, m, n), (q, k, o) = coeffs(*A, z, l)
    for value in (
        p * x * x + m * x + n,
        q * y * y + k * y + o,
        2 * p * x + m + 2 * q * y + k,
        co[0] + co[1] * x - f,
        co[4] + co[5] * y - e,
    ):
        assert sp.cancel(value) == 0
    _, opposite, _ = gram_twist(A, B, z, l)
    for color in range(4):
        indices = [4 * i + (i ^ color) for i in range(4)]
        det = sp.factor(opposite.extract(indices, indices).det())
        assert sp.factor(det + (l - 1) ** 2 * (l * z - 1) ** 2 / l**2) == 0
    field = rational_function_field(("r", "l"))
    sample = {r: 2, l: 3}
    windows = relation_bases(A, B, z, l, sample, field)
    for start in range(4):
        assert PathRow(windows, start, 5, field.one).dims == [1, 4, 10, 20, 36, 68]
    print(sign, "all four generic starts have dimensions 36 and 68", flush=True)


def check_exhaustion(r, l):
    # Away from the poles, det(B)=-det(A) gives f = D/(r*(l^2*z-1))
    # for one of the two choices of r. The remaining curve equation has two roots.
    t = sp.Symbol("t")
    z = r * r
    D = (l * l * z + 1) * t - 2 * l * z
    f = D / (r * (l * l * z - 1))
    co, _, _ = coeffs(t, sp.S.One, z, t, z, l)
    residual = sp.factor(
        z * co[1] * f * f + f - co[0] - co[1] * t * t
    ).as_numer_denom()[0]
    product = D
    specialized = []
    for sign in (1, -1):
        A, _ = section(r, l, sign)
        product *= A[0, 1] * t - A[0, 0]
        pol = sp.resultant((A[0, 1] * t - A[0, 0]).subs(l, 2), r * r - 3, r)
        specialized.append(sp.Poly(pol, t).monic())
    assert not sp.cancel(residual / product).has(t)
    assert sp.cancel(residual / product) != 0
    assert set(specialized) == {
        sp.Poly(t * t - 78 * t + 69, t, domain=sp.QQ),
        sp.Poly(t * t - sp.Rational(78, 23) * t + sp.Rational(3, 23), t, domain=sp.QQ),
    }


if __name__ == "__main__":
    r, l = sp.symbols("r l")
    check_exhaustion(r, l)
    for sign in (1, -1):
        verify_section(sign, r, l)
