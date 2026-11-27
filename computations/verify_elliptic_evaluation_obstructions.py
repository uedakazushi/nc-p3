"""Exact growth obstructions on the full-thin chart with invertible cross blocks."""
from itertools import product

import sympy as sp

from flint_rational_function import rational_function_field
from verify_klein_center_lift import echelon
from verify_klein_elliptic_completion import gram_twist
from verify_parallel_invertible_rectangles import matrix


def relations(a, b, z, scale):
    tensor, _, _ = gram_twist(a, b, z, scale)
    coefficients = {(i, j, k, i ^ j ^ k): tensor[4*i+j, 4*k+(i ^ j ^ k)]
                    for i, j, k in product(range(4), repeat=3)}
    result = []
    for start in range(4):
        edge = []
        for color in range(4):
            cut = matrix(coefficients, start, color)
            for column in range(2 if color < 2 else 1):
                edge.append(sp.Matrix(4, 4, lambda i, j:
                            cut[i, column] if i ^ j == color else 0))
        result.append(edge)
    return result


def kernel_plane(evaluation):
    thin = evaluation[-2:, :].applyfunc(sp.factor)
    plane = (-thin[:, :2].inv()*thin[:, 2:]).col_join(sp.eye(2)).applyfunc(sp.factor)
    assert (evaluation*plane).applyfunc(sp.factor) == sp.zeros(6, 2)
    return plane


def restricted_block(edges, point):
    left = kernel_plane(sp.Matrix.vstack(*[point.T*r for r in edges[0]]))
    right = kernel_plane(sp.Matrix.vstack(*[(r*point).T for r in edges[3]]))
    def flatten(r):
        return sp.Matrix(list(r.applyfunc(sp.factor)))
    unit = sp.eye(2)
    columns = [sp.kronecker_product(flatten(left.T*r), unit[:, i])
               for r in edges[1] for i in range(2)]
    columns += [sp.kronecker_product(unit[:, i], flatten(r*right))
                for r in edges[2] for i in range(2)]
    return sp.Matrix.hstack(*columns)


def exact_rank(mat, parameters):
    field = rational_function_field(tuple(map(str, parameters)))
    rows = [{i: field.from_sympy(x) for i, x in enumerate(row) if x != 0}
            for row in mat.tolist()]
    return len(echelon(rows))


def full_support(s, m, z, scale):
    graph = sp.Matrix([[1, z*m], [m, 1]]).inv()*sp.Matrix([[1, s], [s, 1]])
    d0 = (sp.diag(1, s)*graph*sp.diag(1, 1/m)).T
    d1 = (sp.Matrix([[0, 1], [s, 0]])*graph*sp.Matrix([[0, 1], [m, 0]]).inv()).T
    def cross(l):
        return sp.Matrix([[1, l], [l*z, 1]]), sp.Matrix([[1, 1/l], [1/l, z]])
    y, x = cross(scale)
    a = (d0.inv()*y).applyfunc(sp.factor)
    b = (d1.inv()*scale*x.T).T.applyfunc(sp.factor)
    c = sp.Matrix([[a[0], b[1]], [b[2], a[3]]])
    d = sp.Matrix([[b[0], a[1]], [a[2], b[3]]])
    def schur(corner, l):
        y, x = cross(l)
        return (y*corner.inv()*x).applyfunc(sp.factor)
    u, v = schur(a, scale), schur(b, scale)
    assert (schur(c, sp.S.One)-sp.Matrix([[u[0], v[1]], [v[2], u[3]]])).applyfunc(sp.factor).is_zero_matrix
    assert (schur(d, sp.S.One)-sp.Matrix([[v[0], u[1]], [u[2], v[3]]])).applyfunc(sp.factor).is_zero_matrix
    return a, b


if __name__ == '__main__':
    s, m, z, scale = sp.symbols('s m z l')
    c, d, f, h = sp.symbols('c d f h')
    a, b = sp.Matrix([[1, scale], [c, d]]), sp.Matrix([[scale, f], [1, h]])
    block = restricted_block(relations(a, b, z, scale), sp.Matrix([1, 0, 1, 0]))
    assert exact_rank(block, (c, d, f, h, z, scale)) == 14
    a, b = full_support(s, m, z, scale)
    block = restricted_block(relations(a, b, z, scale), sp.Matrix([1, s, 1, m]))
    assert exact_rank(block, (s, m, z, scale)) == 12
    print('Both Schur completions agree. Three-letter relation ranks are 14 and 12.')
