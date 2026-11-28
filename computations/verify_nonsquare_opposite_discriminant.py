"""Exact polynomial identities for the branched-cover and ternary-block arguments."""

from itertools import combinations_with_replacement

import sympy as sp


def symmetric_product(a, b):
    return (a * b.T + b * a.T) / 2


def equal(left, right):
    return all(sp.expand(entry) == 0 for entry in left - right)


def verify_covers():
    u, w, v, q = sp.symbols("u w v q")
    x, y = u * v, w * q
    z, t = (u * q + w * v) / 2, (u * q - w * v) / 2
    a = sp.Matrix(sp.symbols("a:3"))
    b = sp.Matrix(sp.symbols("b:3"))
    assert sp.expand(t**2 - z**2 + x * y) == 0
    assert equal(symmetric_product(u * a + w * b, v * a + q * b),
                 x * a * a.T + 2 * z * symmetric_product(a, b) + y * b * b.T)
    assert equal((u * a + w * b).cross(v * a + q * b), 2 * t * a.cross(b))
    assert equal(symmetric_product(u * a + w * b, u * a - w * b),
                 u**2 * a * a.T - w**2 * b * b.T)
    second = symmetric_product(sp.Matrix([u**2, u*w, w**2]), sp.Matrix([v**2, v*q, q**2]))
    expected = sp.Matrix([[x**2, x*z, 2*z**2-x*y], [x*z, x*y, y*z],
                          [2*z**2-x*y, y*z, y**2]])
    assert equal(second, expected)


def verify_ternary_block():
    d = sp.Matrix(sp.symbols("d:3"))
    i, j, k = d
    cross = sp.Matrix([[0, -k, j], [k, 0, -i], [-j, i, 0]])
    c = sp.diag(*sp.symbols("c:3"))
    h = cross * c * cross.T
    assert equal(h.adjugate(), (d.T * c.adjugate() * d)[0] * d * d.T)
    h = cross * cross.T
    pairs = list(combinations_with_replacement(range(3), 2))
    coefficients = [[sp.Poly(h[i, j], *d).coeff_monomial(d[a] * d[b])
                     for a, b in pairs] for i, j in pairs]
    assert sp.Matrix(coefficients).rank() == 6
    v = sp.Matrix(sp.symbols("v:4"))
    a, b, c, f, g, h = sp.symbols("a b c f g h")
    block = sp.Matrix([[a, b, c], [b, f, g], [c, g, h]])
    matrix = v * v.T + sp.diag(block, 0)
    assert equal(matrix.adjugate()[:3, :3], v[3]**2 * block.adjugate())


if __name__ == "__main__":
    verify_covers()
    verify_ternary_block()
    print("Verified the two cover models, the contact quadrics, and both adjugate formulas.")
