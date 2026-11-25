"""Characteristic-zero ranks on both entire smooth Hesse torsion curves."""

from itertools import product

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from ternary_ore_parameters import parameters


def block(a, b, cs):
    # Weighted diagonal (1,omega,omega^2): degree-zero summand of Hom(H,B2).
    base = range(3)
    eye = sp.eye(9)
    relations = sp.Matrix.hstack(
        *[
            cs[j] * eye[:, 4 * j]
            + a * eye[:, 3 * ((j + 1) % 3) + (j + 2) % 3]
            + b * eye[:, 3 * ((j + 2) % 3) + (j + 1) % 3]
            for j in base
        ]
    )
    pure = sp.kronecker_product(relations, sp.eye(3)).row_join(
        sp.kronecker_product(sp.eye(3), relations)
    )
    lifts = []
    for j in base:
        for i, k in product(base, repeat=2):
            if i != k and (i + k - j) % 3 == 0:
                op = sp.zeros(9, 3)
                op[3 * i + k, j] = 1
                lifts.append(op)
    left, right = [], []
    denominators = []
    for j in base:
        rows = [
            9 * l + 3 * m + n
            for l, m, n in product(base, repeat=3)
            if (l + m + n - 2 * j) % 3 == 0
        ]
        for normal in pure[rows, :].T.nullspace():
            normal = normal.applyfunc(sp.factor)
            denominators.extend(v.as_numer_denom()[1] for v in normal)
            for into, first in ((left, True), (right, False)):
                row = []
                for op in lifts:
                    action = (
                        sp.kronecker_product(op, sp.eye(3))
                        if first
                        else sp.kronecker_product(sp.eye(3), op)
                    )
                    row.append(
                        sp.factor((normal.T * (action * relations[:, j])[rows, :])[0])
                    )
                into.append(row)
    return (
        sp.Matrix(left),
        sp.Matrix(right),
        sp.factor(sp.lcm_list(denominators)),
        relations,
    )


def main():
    t = sp.symbols("t", nonzero=True)
    for n, b, u in [(2, 1, t), (4, t, t * (1 + t**2) / (1 + t))]:
        left, right, den, rels = block(1, b, (1, 1, u))
        assert left.shape == right.shape == (10, 6)
        cols = [0, 1, 2, 3, 5, 6]
        rows = [0, 1, 2, 3, 4, 6] if n == 2 else list(range(6))
        det = sp.factor(
            left.row_join(right).extract(rows, cols).det(method="domain-ge")
        )
        expected = (
            -t * (t - 1)
            if n == 2
            else -(t**2) * (t - 1) * (t**2 + 1) * (t**2 + t + 1) / (t + 1) ** 2
        )
        assert sp.cancel(det - expected) == 0
        # Only poles of the actual quotient functionals (including their normal vectors)
        # are allowed in the uniform lower bound from this minor.
        excluded = sp.factor(t * u * (1 + b) * ((1 + b**3 + u) ** 3 - 27 * b**3 * u))
        bad = sp.factor(
            excluded.as_numer_denom()[0] * excluded.as_numer_denom()[1] * (t + 1)
        )
        for expr in [den, det.as_numer_denom()[0]]:
            for f, _ in sp.factor_list(expr)[1]:
                assert sp.rem(bad, f, t) == 0
        eigenspaces = {1: 4, -1: 3} if n == 2 else {-1: 4, sp.I: 4, -sp.I: 4}
        for mu, rnk in eigenspaces.items():
            assert DomainMatrix.from_Matrix(left + mu * right).rank() == rnk
        dims, gauge, shift = parameters(rels.subs(t, 2))
        assert dims == [18, 15] and gauge == 6 and shift**n == sp.eye(9)
        print(
            f"order {n}, minor {det}, denominator {den}, generic pencil ranks {eigenspaces}"
        )
    print(
        "PASS: two full smooth Hesse loci and the nine-dimensional shifts at rational fibers"
    )


if __name__ == "__main__":
    main()
