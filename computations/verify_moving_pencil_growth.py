"""Check double-diagonal path evaluation, including coincident successors."""
from sympy import Matrix


def multiply_increment(poly, variable):
    """Multiply by t_0 + ... + t_variable with every t_i squared zero."""
    result = {}
    for monomial, coefficient in poly.items():
        for index in range(variable + 1):
            bit = 1 << index
            if monomial & bit:
                continue
            term = monomial | bit
            result[term] = result.get(term, 0) + coefficient
    return result


def verify_path_evaluation(n):
    # Fix x_0=0 and write x_j as the sum of the first j increments.
    columns = [{0: 1}]
    for variable in range(n):
        columns += [multiply_increment(p, variable) for p in columns]
    # In subset order the product for S has leading mask S with coefficient 1.
    for mask, column in enumerate(columns):
        assert max(column) == mask and column[mask] == 1
    if n <= 5:
        size = 1 << n
        assert Matrix(size, size, lambda i,j: columns[j].get(i, 0)).det() == 1


for length in range(1, 9):
    verify_path_evaluation(length)

# Ranks and determinant degrees on P^1 x P^1, for each projection.
def external(left, right):
    rank, degree = left
    other_rank, other_degree = right
    return rank * other_rank, degree * other_rank, rank * other_degree


quotient, relations, middle = (2, 1), (4, 1), 4
target = tuple(middle * x for x in external(quotient, quotient))
domain = tuple(a+b for a,b in zip(external(relations, quotient), external(quotient, relations)))
assert target == (16, 8, 8) and domain == (16, 6, 6)
print("PASS: determinant bidegree (2,2) and independent binary words on double paths through length 8")
