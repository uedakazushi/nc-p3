"""Check all mixed-exchange fibers, including the zero-trace and eight exceptional charts."""
from itertools import product

import sympy as sp

from verify_elliptic_evaluation_obstructions import exact_rank, restricted_block
from verify_klein_elliptic_completion import gram_twist
from verify_parallel_invertible_rectangles import matrix


def edges(a, b, d, f, g, scale):
    first, second = sp.Matrix([[a, b], [b, d]]), sp.Matrix([[b, f], [g, b]])
    tensor, _, _ = gram_twist(first, second, sp.S.One, scale)
    coefficients = {(i, j, k, i ^ j ^ k): tensor[4*i+j, 4*k+(i ^ j ^ k)]
                    for i, j, k in product(range(4), repeat=3)}
    result = []
    for start in range(4):
        edge = []
        for color in range(4):
            cut = matrix(coefficients, start, color)
            columns = (0,) if color >= 2 else ((2, 3) if start%2 and color == 1 else (0, 1))
            for column in columns:
                edge.append(sp.Matrix(4, 4, lambda i, j:
                            sp.factor(cut[i, column]) if i ^ j == color else 0))
        result.append(edge)
    return result


def units(expressions):
    result = set()
    for expression in expressions:
        numerator, denominator = sp.fraction(sp.factor(expression))
        assert numerator != 0
        for polynomial in (numerator, denominator):
            result.update(factor for factor, _ in sp.factor_list(polynomial)[1])
    return result


def good_chart(trace, b, scale, u, v):
    k = (trace*b-scale)*(trace*b-1/scale)
    if trace:
        a, d, f, g = 1+u, 1+k/u, 1+v, 1+k/v
    else:
        a, d, f, g = u, (scale+1/scale)*b-u, v, (scale+1/scale)*b-v
    allowed = units((scale, scale**2-1, k, a*d-f*g, a*d-b*b, b*b-f*g))
    relations = edges(a, b, d, f, g, scale)
    for start in (1, 3):
        for color in (0, 1):
            frame = sp.Matrix(4, 2, lambda i, j: relations[start][2*color+j][i, i ^ color])
            corner = frame[2:, :] if color else frame[:2, :]
            other = frame[:2, :] if color else frame[2:, :]
            assert other == sp.ones(2)
            assert units([corner.det(), *(other[0, :]*corner.inv())]) <= allowed
    y = sp.Matrix([[1, scale], [scale, 1]])
    assert units(list(y*sp.Matrix([[a, b], [b, d]]).inv())) <= allowed


def exceptional(scale, u, v):
    for i, j, inv in product(range(2), repeat=3):
        a, d = (u, 1) if i else (1, u)
        f, g = (v, 1) if j else (1, v)
        b = 1/scale if inv else scale
        rel = edges(a, b, d, f, g, scale)
        start = (1+inv) if i == j else (3 if inv else 0)
        point = sp.eye(4)[:, i]+sp.eye(4)[:, i+2]
        block = restricted_block(rel[start:]+rel[:start], point)
        assert exact_rank(block, (scale, u, v)) == 12


if __name__ == '__main__':
    b, scale, u, v = sp.symbols('b l u v')
    for trace in (0, 1):
        good_chart(trace, b, scale, u, v)
    exceptional(scale, u, v)
    print('Both regular charts have full graph weights. All eight exceptional block ranks are 12.')
