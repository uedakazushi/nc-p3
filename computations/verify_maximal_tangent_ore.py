"""A finite symmetric AS web and the cubic Ore conditions at its diagonal point."""

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from ternary_ore_parameters import parameters


def rank(matrix):
    return DomainMatrix.from_Matrix(matrix).rank()


def columns(vectors):
    return sp.Matrix.hstack(*vectors)


def overlap(left, right, n):
    return sp.kronecker_product(left, sp.eye(n)).row_join(
        sp.kronecker_product(sp.eye(n), right)
    )


def main():
    x, y, z, w = variables = sp.symbols("x y z w")
    web = [x * x + 4 * y * z + 2 * w * x, y * y + 4 * x * z, z * z + 4 * x * y, w * w]
    matrices = [sp.hessian(q, variables) / 2 for q in web]
    assert sp.groebner(web, *variables).is_zero_dimensional
    s, t, u = params = sp.symbols("s t u")
    h = sum((p * m[:3, :3] for p, m in zip(params, matrices[:3])), sp.zeros(3))
    v = sp.Matrix([s, 0, 0])
    det = h.det()
    assert sp.groebner([det.diff(p) for p in params], *params).is_zero_dimensional
    # No rank-one ternary member; away from w^2 the upper block must have rank two.
    # For rank two the lift exists precisely on the zero locus of adj(h)*v.
    constraints = [det, *(h.adjugate() * v)]
    for i, p in enumerate(params):
        subs = {q: 0 for q in params[:i]} | {p: 1}
        remaining = params[i + 1 :]
        if remaining:
            gb = sp.groebner([f.subs(subs) for f in constraints], *remaining)
            assert list(gb) == [1] or gb.is_zero_dimensional
        else:
            assert any(f.subs(subs) for f in constraints)
    pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    basis = []
    for i, j in pairs:
        b = sp.zeros(4)
        b[i, j] = b[j, i] = 1
        basis.append(b)
    pairing = sp.Matrix(
        [[sum(b.multiply_elementwise(m)) for b in basis] for m in matrices]
    )
    relations = [
        sum((a * b for a, b in zip(vec, basis)), sp.zeros(4))
        for vec in pairing.nullspace()
    ]
    r = columns([q.reshape(16, 1) for q in relations])
    contraction = columns([q[3, :3].T for q in relations])
    assert rank(contraction) == 3 and all(q[3, 3] == 0 for q in relations)
    jacobian = columns([q[:3, 3].col_join(q[3, :3].T) for q in relations])
    assert rank(jacobian) == 3
    pure = r * columns(contraction.nullspace())
    hrows = [4 * i + j for i in range(3) for j in range(3)]
    ternary = pure[hrows, :]
    assert pure[[i for i in range(16) if i not in hrows], :] == sp.zeros(7, 3)
    cubic = overlap(ternary, ternary, 3)
    assert rank(cubic) == 17 and 64 - rank(overlap(r, r, 4)) == 20
    syzygy = cubic[:, :9] * cubic.nullspace()[0][:9, :]
    assert rank(syzygy.reshape(3, 9)) == rank(syzygy.reshape(9, 3)) == 3
    dimensions, translations, shift = parameters(ternary)
    assert dimensions == [18, 15] and translations == 6
    assert shift**2 == sp.eye(9) and rank(shift - sp.eye(9)) == 6
    print(
        "AS web: core dimensions",
        dimensions,
        "translations",
        translations,
        "shift characteristic polynomial",
        shift.charpoly().as_expr().factor(),
    )
    for sign in (-1, 1):
        basis = sp.eye(9)
        standard = columns(
            [
                basis[:, 3 * i + j] + sign * basis[:, 3 * j + i]
                for i in range(3)
                for j in range(i + 1, 3)
            ]
        )
        dims, gauges, action = parameters(standard)
        if sign == -1:
            assert dims == [18] and gauges == 3 and action == sp.eye(15)
        else:
            assert dims == [18, 15] and gauges == 6
            assert action**2 == sp.eye(9) and rank(action - sp.eye(9)) == 6
        print(
            "base sign",
            sign,
            "core dimensions",
            dims,
            "translations",
            gauges,
            "shift characteristic polynomial",
            action.charpoly().as_expr().factor(),
        )
    # Changing one mixed correction at one start leaves its tangent and its
    # pure ternary relations unchanged but can destroy cubic compatibility.
    mixed = next(i for i in range(r.cols) if contraction[:, i] != sp.zeros(3, 1))
    # Change by a multiple of the contraction in every relation, fixing the pure kernel.
    a = next(j for j in range(3) if contraction[j, mixed])
    altered = r.copy()
    altered[0, :] += contraction[a, :]
    assert 64 - rank(overlap(r, altered, 4)) < 20
    print(
        "PASS: finite AS web, three-dimensional tangent, unique ternary syzygy of end ranks 3, cubic failure after a mixed perturbation"
    )


if __name__ == "__main__":
    main()
