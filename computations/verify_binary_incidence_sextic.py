"""A smooth residual plane sextic for a lifted binary relation space."""

import sympy as sp


def convolution(degree, coeffs):
    return sp.Matrix(
        degree + 1, degree + 3, lambda i, j: coeffs[j - i] if 0 <= j - i <= 2 else 0
    )


def check():
    coords = sp.symbols("a b c")
    lift = sp.Matrix(
        [
            [1, -7, -2, -6, -7, -4],
            [-5, 6, 6, 3, 8, -8],
            [5, 8, -6, 8, 0, -2],
            [-7, 4, -8, -4, 2, 0],
        ]
    )
    matrix = convolution(3, coords).col_join(convolution(1, coords) * lift)
    form = matrix.det(method="domain-ge")
    poly = sp.Poly(form, *coords)
    assert poly.is_homogeneous and poly.total_degree() == 6
    gradients = [sp.diff(form, q) for q in coords]
    for i in range(3):
        patch = {q: 0 for q in coords[:i]}
        patch[coords[i]] = 1
        eqs = [f.subs(patch) for f in gradients]
        if i < 2:
            assert list(
                sp.groebner(eqs, *coords[i + 1 :], modulus=1009, order="grevlex")
            ) == [1]
        else:
            assert any(int(c) % 1009 for c in eqs)
    print(
        "The plane sextic determined by the displayed integral lift is smooth modulo 1009 and in characteristic zero."
    )


if __name__ == "__main__":
    check()
