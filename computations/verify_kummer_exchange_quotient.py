"""Exclude fixed points of the resolved block exchange on one proper fiber."""

from itertools import combinations

import sympy as sp

from verify_kummer_discriminant import discriminant_form
from verify_kummer_polar import polarity


def check():
    variables, form = discriminant_form()
    z, l = variables[-2:]
    coords = variables[:4]
    form = sp.expand(form.subs({z: 3, l: 2}))
    matrix = polarity(sp.Integer(3), sp.Integer(2))
    prime = 1009
    assert matrix.det() % prime
    x = sp.Matrix(coords)
    polar = matrix * sp.Matrix([sp.diff(form, u) for u in coords])
    isotropic = (x.T * matrix.adjugate() * x)[0]
    minors = [polar[i] * x[j] - polar[j] * x[i] for i, j in combinations(range(4), 2)]
    for i in range(4):
        patch = {u: 0 for u in coords[:i]}
        patch[coords[i]] = 1
        generators = [sp.expand(f.subs(patch)) for f in [form, isotropic, *minors]]
        if i == 3:
            assert any(int(f) % prime for f in generators)
        else:
            basis = sp.groebner(
                generators, *coords[i + 1 :], modulus=prime, order="grevlex"
            )
            assert list(basis) == [1]
    print(
        "The projective fixed-point candidate scheme is empty at (3,2) modulo 1009, including singular points."
    )


if __name__ == "__main__":
    check()
