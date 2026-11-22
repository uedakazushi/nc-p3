"""All three tangent defects on the smooth curve y^2=x^3-x over Q."""
from itertools import combinations, product

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



def marked_jacobian(beta, metric):
    """Differentiate the mixed Gram equations, fixing both marked axes."""
    rows = []
    pairs = list(combinations(range(4), 2))
    for (i, k), (j, l) in product(pairs, repeat=2):
        if i == 0 and j == 0:
            continue
        row = sp.zeros(1, 54)
        for a, b, c, d, sign in [(i,j,k,l,1),(k,l,i,j,1),
                                  (i,l,k,j,-1),(k,j,i,l,-1)]:
            if a and b:
                start = 6*(3*(a-1)+b-1)
                row[:, start:start+6] += sign*beta[c,d].T*metric
        rows.append(row)
    return sp.Matrix.vstack(*rows)


def marked_control(first, second, left, right, middle, pairing):
    # Lift quotient bases to H^0(N), then discard the two pencil coordinates.
    pool = [first, second] + middle
    basis = [pool[i] for i in coordinates(pool).rref()[1]]
    domain = [first]+left
    basis_left = [domain[i] for i in coordinates(domain).rref()[1]]
    values = [a*b for a, b in product(basis_left, right)]
    data = coordinates(basis+values)
    lifts, _ = data[:, :8].gauss_jordan_solve(data[:, 8:])
    beta = {(i,j): lifts[2:,4*i+j] for i,j in product(range(4), repeat=2)}
    metric = pairing(basis[2:], basis[2:])
    assert metric.rank() == 6
    jet = marked_jacobian(beta, metric)
    u = sp.Matrix.hstack(*(beta[0,j] for j in (1,2,3)))
    radical = (u.T*metric*u).nullspace()
    ranks = []
    for vector in radical:
        r = u*vector
        c = sp.Matrix(3, 3, lambda i,j: (r.T*metric*beta[i+1,j+1])[0])
        ranks.append(c.rank())
    return jet.rank(), ranks



def constant_rank_pencils():
    s, t = sp.symbols("s t")
    common_kernel = sp.Matrix([[s, 0, 0], [t, s, 0], [0, t, 0]])
    moving_kernels = sp.Matrix([[0, 0, s], [0, 0, t], [s, t, 0]])
    pairs = list(combinations(range(3), 2))
    for pencil, expected in [(common_kernel, 9), (common_kernel.T, 9), (moving_kernels, 8)]:
        columns = []
        for coefficient in [pencil.diff(s), pencil.diff(t)]:
            for i, j in product(range(3), repeat=2):
                h = sp.zeros(3)
                h[i,j] = 1
                values = [coefficient[a,b]*h[k,l] + h[a,b]*coefficient[k,l]
                          - coefficient[a,l]*h[k,b] - h[a,l]*coefficient[k,b]
                          for (a,k), (b,l) in product(pairs, repeat=2)]
                columns.append(sp.Matrix(values))
        jacobian = sp.Matrix.hstack(*columns)
        assert jacobian.rank() == expected
        if expected == 8:
            assert jacobian.T.nullspace() == [sp.eye(9)[:, 0]]
    print("constant-rank radical pencils: cofactor differential ranks 9, 9, 8")


def main():
    constant_rank_pencils()
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
        jet_rank, radical_ranks = marked_control(first, second, left, right, middle, pairing)
        assert jet_rank == 27
        if defect == 1:
            assert radical_ranks == [2]
        print("marked section", first, "defect", defect, "Gram rank", q.rank(),
              "mixed Jacobian rank", jet_rank, "radical form ranks", radical_ranks)


if __name__ == "__main__":
    main()
