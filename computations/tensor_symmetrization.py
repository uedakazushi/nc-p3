"""Exact linear symmetrization and rational invertible-witness construction."""
from fractions import Fraction as F
from sympy import Matrix, symbols, expand, zeros
from frobenius_tensor import digits, nullspace_sparse


def encode(word):
    result = 0
    for value in word:
        result = 4 * result + value
    return result


def act_factor(tensor, matrix, factor):
    result = {}
    for key, value in tensor.items():
        word = list(digits(key))
        source = word[factor]
        for target in range(4):
            coefficient = F(matrix[target, source])
            if coefficient:
                word[factor] = target
                index = encode(word)
                result[index] = result.get(index, F(0)) + coefficient * value
    return {i: x for i, x in result.items() if x}


def swap_factors(tensor, left, right):
    result = {}
    for key, value in tensor.items():
        word = list(digits(key))
        word[left], word[right] = word[right], word[left]
        result[encode(word)] = value
    return result


def symmetrizer_space(tensor, factor):
    swapped = swap_factors(tensor, 0, factor)
    equations = {}
    for column in range(16):
        elementary = zeros(4)
        elementary[column // 4, column % 4] = 1
        lhs = act_factor(tensor, elementary, factor)
        rhs = act_factor(swapped, elementary, 0)
        for key in lhs.keys() | rhs.keys():
            value = lhs.get(key, F(0)) - rhs.get(key, F(0))
            if value:
                equations.setdefault(key, {})[column] = value
    kernel = nullspace_sparse(equations.values(), 16)
    return [Matrix(4, 4, lambda i, j: vector.get(4 * i + j, F(0)))
            for vector in kernel]


def invertible_witness(basis):
    if not basis:
        return None
    for matrix in basis:
        if matrix.det():
            return matrix
    variables = symbols(f't0:{len(basis)}')
    family = sum((t * matrix for t, matrix in zip(variables, basis)), zeros(4))
    determinant = expand(family.det(method='domain-ge'))
    if determinant == 0:
        return None
    chosen = {}
    for variable in variables:
        for value in range(5):
            remaining = expand(determinant.subs(variable, value))
            if remaining != 0:
                chosen[variable] = value
                determinant = remaining
                break
        else:
            raise AssertionError('nonzero degree-four determinant vanished on grid')
    result = family.subs(chosen)
    assert result.det() != 0
    return result


def symmetrize(tensor):
    spaces = [symmetrizer_space(tensor, j) for j in (1, 2, 3)]
    witnesses = [invertible_witness(space) for space in spaces]
    if any(matrix is None for matrix in witnesses):
        return None, spaces, witnesses
    result = tensor
    for factor, matrix in enumerate(witnesses, 1):
        result = act_factor(result, matrix, factor)
    assert all(swap_factors(result, 0, j) == result for j in (1, 2, 3))
    return result, spaces, witnesses
