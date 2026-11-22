#!/usr/bin/env python3
"""Local equations of the reversible rank-six tensor at the reduced symmetric example."""
from itertools import combinations, combinations_with_replacement

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from verify_finite_window_rigidity import verify_flip_obstructions


def alternation(tensor):
    pairs = list(combinations(range(4), 2))
    return sp.Matrix([tensor[4*a+i, 4*b+j] - tensor[4*a+j, 4*b+i]
                      for a, b in pairs for i, j in pairs])


def opposite(tensor):
    pairs = list(combinations_with_replacement(range(4), 2))
    return sp.Matrix([[tensor[4*a+i, 4*b+j] for i, j in pairs] for a, b in pairs])


def local_matrices():
    rel, h = verify_flip_obstructions()
    columns = []
    for j in range(6):
        for vector in sp.eye(16).columnspace():
            term = rel[:, j]*vector.T
            columns.append(alternation(term + term.T))
    derivative = sp.Matrix.hstack(*columns)
    dependency, = derivative.T.nullspace()
    tangent = sp.Matrix.hstack(*derivative.nullspace())
    assert tangent.cols == 61

    form = sp.zeros(16)
    pairs = list(combinations(range(4), 2))
    for index, ((a, b), (i, j)) in enumerate((p, q) for p in pairs for q in pairs):
        value = dependency[index]/2
        form[4*a+i, 4*b+j] += value
        form[4*b+j, 4*a+i] += value
        form[4*a+j, 4*b+i] -= value
        form[4*b+i, 4*a+j] -= value
    restricted = tangent.T*sp.kronecker_product(h.inv(), form)*tangent
    assert DomainMatrix.from_Matrix(restricted).rank() == 6
    base = opposite(rel*h*rel.T)
    left = sp.Matrix.hstack(*base.T.nullspace())
    right = sp.Matrix.hstack(*base.nullspace())
    assert right.cols == 4

    for index in range(tangent.cols):
        if restricted[index, index] != 0 or restricted[:, index].is_zero_matrix:
            continue
        direction = tangent[:, index].reshape(6, 16).T
        first = opposite(rel*direction.T + direction*rel.T)
        if (left.T*first*right).rank() >= 2:
            return rel, h, derivative, dependency, tangent, restricted, index
    raise AssertionError("No smooth isotropic direction with opposite rank increment two")


def verify_opposite_normal(rel, h, tangent, quadratic):
    base = opposite(rel*h*rel.T)
    assert base == base.T
    kernel = sp.Matrix.hstack(*base.nullspace())
    pairs = list(combinations(range(4), 2))
    columns = []
    for i in range(tangent.cols):
        y = tangent[:, i].reshape(6, 16).T
        block = kernel.T*opposite(rel*y.T + y*rel.T)*kernel
        assert block == -block.T
        columns.append(sp.Matrix([block[a, b] for a, b in pairs]))
    linear = sp.Matrix.hstack(*columns)
    pfaffian = sp.zeros(6)
    for a, b, sign in ((0, 5, 1), (1, 4, -1), (2, 3, 1)):
        pfaffian[a, b] = pfaffian[b, a] = sp.Rational(sign, 2)
    assert linear.rank() == 6
    assert linear.T*pfaffian*linear == sp.Rational(15, 16)*quadratic

    variables = sp.symbols("z:6")
    skew = sp.zeros(4)
    for (a, b), v in zip(pairs, variables):
        skew[a, b], skew[b, a] = v, -v
    minors = [sp.Poly(skew.extract(row, col).det(), *variables) for row in pairs for col in pairs]
    monomials = [a*b for a, b in combinations_with_replacement(variables, 2)]
    assert sp.Matrix([[minor.coeff_monomial(m) for m in monomials] for minor in minors]).rank() == 21


def verify_symmetries(rel, h, tangent):
    # Enumerate all planes containing seven points and check the symmetrizer of R.
    from verify_finite_square_chart import frame
    points, _ = frame()
    large_planes = set()
    for triple in combinations(range(20), 3):
        equations = points[list(triple), :].nullspace()
        if len(equations) != 1:
            continue
        zeros = tuple(i for i, v in enumerate(points*equations[0]) if v == 0)
        if len(zeros) >= 7:
            large_planes.add(zeros)
    assert len(large_planes) == 10 and all(len(p) == 7 for p in large_planes)

    units = sp.eye(16).columnspace()
    symmetrizers = []
    for unit in units:
        end = unit.reshape(4, 4)
        constraints = []
        for r in rel.columnspace():
            matrix = r.reshape(4, 4)*end.T
            constraints.extend(matrix[a, b] - matrix[b, a] for a, b in combinations(range(4), 2))
        symmetrizers.append(sp.Matrix(constraints))
    assert sp.Matrix.hstack(*symmetrizers).rank() == 15

    base = rel*h*rel.T
    orbit = []
    for unit in units:
        end = unit.reshape(4, 4)
        for action in (sp.kronecker_product(end, sp.eye(4)), sp.kronecker_product(sp.eye(4), end)):
            orbit.append(opposite(action*base+base*action.T).reshape(100, 1))
    orbit = sp.Matrix.hstack(*orbit)
    def modular_rank(matrix):
        reduced = matrix.applyfunc(lambda v: int(v.p)*pow(int(v.q), -1, 101) % 101)
        return DomainMatrix.from_Matrix(reduced).convert_to(sp.GF(101)).rank()
    assert modular_rank(orbit) == 31

    frame_points = [*sp.eye(4).columnspace(), -sp.ones(4, 1)]
    exchange = sp.Matrix(16, 16, lambda i, j: int(j == 4*(i % 4) + i//4))
    # Conjugacy classes of the diagonal S5. For each, its composition with
    # exchange acts nontrivially on the tangent space modulo the change of bases.
    for cycle_type in ((1,1,1,1,1), (2,1,1,1), (2,2,1), (3,1,1), (3,2), (4,1), (5,)):
        perm = list(range(5))
        start = 0
        for size in cycle_type:
            for offset in range(size):
                perm[start+offset] = start+(offset+1) % size
            start += size
        action = sp.Matrix.hstack(*(frame_points[i] for i in perm[:4])).inv().T
        swap = sp.kronecker_product(action, action)*exchange
        assert swap*base*swap.T == base
        for column in range(tangent.cols):
            y = tangent[:, column].reshape(6, 16).T
            variation = y*rel.T + rel*y.T
            difference = opposite(swap*variation*swap.T-variation).reshape(100, 1)
            if modular_rank(orbit.row_join(difference)) == 32:
                break
        else:
            raise AssertionError(cycle_type)


if __name__ == "__main__":
    rel, h, _, _, tangent, quadratic, _ = local_matrices()
    verify_opposite_normal(rel, h, tangent, quadratic)
    verify_symmetries(rel, h, tangent)
    print("PASS derivative rank 35, normal quadratic rank 6, opposite rank at least 8,")
    print("symmetrizer rank 15, ten seven-point planes, all seven period-one tangent exclusions, and the skew normal block with its 21 minor quadrics")
