"""All three tangent defects on the smooth curve y^2=x^3-x over Q."""
from itertools import product

import sympy as sp

x, y = sp.symbols("x y")


def coordinates(expressions):
    """Common Laurent monomial coordinates in Q[x,x^-1]+y Q[x,x^-1]."""
    columns = []
    for expression in expressions:
        terms = {}
        for term in sp.expand(expression).as_ordered_terms():
            coefficient, monomial = term.as_coeff_Mul()
            powers = monomial.as_powers_dict()
            i, j = int(powers[x]), int(powers[y])
            for k in range(j//2+1):
                exponent = (i + j//2 + 2*k, j % 2)
                value = coefficient * sp.binomial(j//2, k) * (-1)**(j//2-k)
                terms[exponent] = terms.get(exponent, 0) + value
        columns.append(terms)
    rows = sorted(set().union(*columns))
    return sp.Matrix([[column.get(row, 0) for column in columns] for row in rows])


def main():
    # L=4O, M=3O+P with P=(0,0); the second marked section is 1 in M.
    left = [1, x, y, x**2]
    right = [1, x, y, y/x]
    mixed = [a*b for a, b in product(left, right)]
    middle = [mixed[i] for i in coordinates(mixed).rref()[1]]
    assert len(middle) == 8
    second = x**2*y + y/x + 11
    # Maximal poles at O and P. Check the other zeroes of the marked sections.
    assert sp.resultant(y**2-6, second.subs(x, 2), y)
    assert sp.resultant(x**3-x-1, x*second.subs(y, -1), x)

    for first, defect in [(1, 2), (x-2, 1), (y+1, 0)]:
        multiples = [a*b for a, b in product([first, second], middle)]
        symmetric = [a*b for a, b in product(middle, repeat=2)]
        all_products = coordinates(multiples+symmetric)
        quotient = all_products[:, :len(multiples)]
        assert quotient.rank() == 15 and all_products.rank() == 16
        trace = next(v for v in quotient.T.nullspace() if (v.T*all_products).rank())

        def pairing(one, two):
            values = [a*b for a, b in product(one, two)]
            full = coordinates(multiples + symmetric + values)
            # With other column sets, globally zero coordinates can appear.
            # Factor through the same rank-16 space using independent columns.
            base = full[:, :all_products.cols]
            indices = all_products.rref()[1]
            coefficients, _ = base[:, indices].gauss_jordan_solve(full[:, -len(values):])
            entries = trace.T*all_products[:, indices]*coefficients
            return entries.reshape(len(one), len(two))

        q = pairing([first*b for b in right], [first*b for b in right])
        jacobian = coordinates([first, second] + left + [first*b for b in right])
        assert 8-jacobian.rank() == 3-q.rank() == defect
        assert pairing(mixed, mixed).rank() == 6
        assert pairing([a*b for a,b in product(left, repeat=2)],
                       [a*b for a,b in product(right, repeat=2)]).rank() == 8
        print("marked section", first, "tangent defect", defect, "Gram rank", q.rank())


if __name__ == "__main__":
    main()
