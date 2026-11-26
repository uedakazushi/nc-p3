"""A full-rank mixed bracket reconstruction in the smooth marked component."""
import numpy as np
import sympy as sp
from verify_mixed_component_smoothing import kernel
from verify_semisimple_period_obstruction import G, P, Q, inverse, reduced

DIGITS = 72
MODULUS = P**DIGITS
INDICES = [(i, j) for i, j in np.ndindex(3, 3) if (i, j) != (2, 2)]


def units(modulus):
    identity = np.eye(3, dtype=object)
    return np.array([identity[:, i:i+1]@identity[j:j+1, :]
                     -int(i == j)*pow(3, -1, modulus)*identity for i, j in INDICES]) % modulus


def unpack(vector, basis):
    return (vector.reshape(4, 3, 8)@basis.reshape(8, 9)).reshape(4, 3, 3, 3)


def brackets(left, right):
    return left[:, None]@right[None, :]-right[None, :]@left[:, None]


def product(left, right):
    return (brackets(left[0], right[1])+brackets(left[2], right[3]).transpose(2, 3, 0, 1)).reshape(-1)


def lift():
    seed = np.array([G, -G, Q, Q], dtype=object)
    small = units(P).astype(np.int64)
    modular_seed = (seed % P).astype(np.int64)
    directions = [unpack(e, small) for e in np.eye(96, dtype=np.int64)]
    jac = np.array([(product(modular_seed, d)+product(d, modular_seed)) % P for d in directions]).T
    rows, cols = reduced(jac.T)[1], reduced(jac)[1]
    assert len(cols) == len(rows) == 62
    correction = inverse(jac[np.ix_(rows, cols)])
    rng = np.random.default_rng(536)
    free = rng.integers(P, size=(2, 34))[1]
    basis = units(MODULUS)
    first = (free@kernel(jac) % P).astype(object)
    current = (seed+P*unpack(first, basis)) % MODULUS
    for n in range(2, DIGITS):
        residual = product(current, current)
        assert not np.any(residual % P**n)
        rhs = (residual[rows]//P**n % P).astype(np.int64)
        delta = np.zeros(96, dtype=object)
        delta[cols] = -correction@rhs % P
        current = (current+P**n*unpack(delta, basis)) % MODULUS
    assert not np.any(product(current, current) % MODULUS)
    return current


def kernel_over_dvr(matrix, precision):
    """Unit-pivot elimination, discarding only digits lost in common divisions."""
    array = np.array(matrix, dtype=object) % P**precision
    n, m = array.shape
    free, pivots = list(range(m)), []
    for row in range(n):
        while not np.any(array[row:, free] % P):
            if not np.any(array[row:, free]):
                raise AssertionError('insufficient precision or dependent rows')
            array[row:, :] //= P
            precision -= 1
            array %= P**precision
        i, j = np.argwhere(array[row:, free] % P != 0)[0]
        column = free.pop(j)
        array[[row, row+i]] = array[[row+i, row]]
        array[row] = array[row]*pow(int(array[row, column]), -1, P**precision) % P**precision
        factors = array[:, column].copy()
        factors[row] = 0
        array = (array-factors[:, None]*array[row]) % P**precision
        pivots.append(column)
    basis = np.zeros((m-n, m), dtype=object)
    basis[:, free] = np.eye(m-n, dtype=object)
    basis[:, pivots] = -array[:, free].T
    return basis % P**precision, precision


def trace_matrix(c, basis):
    traces = (c.reshape(9, 9)@basis.transpose(0, 2, 1).reshape(8, 9).T).reshape(3, 3, 8)
    rows = []
    for j in range(3):
        for i in range(3):
            for k in range(i, 3):
                row = np.zeros((3, 8), dtype=object)
                row[i] += traces[k, j]
                row[k] += traces[i, j]
                rows.append(row.reshape(-1))
    return rows


def check_inverse(c, basis):
    left, n = kernel_over_dvr(trace_matrix(c, basis), DIGITS)
    right, m = kernel_over_dvr(trace_matrix(c.swapaxes(0, 1), basis), DIGITS)
    assert left.shape == right.shape == (6, 24)
    def matrices(rows):
        return (rows.reshape(-1, 3, 8)@basis.reshape(8, 9)).reshape(-1, 3, 3, 3) % P**min(n, m)
    product_rows = np.array([brackets(r, s).reshape(-1) for r in matrices(left) for s in matrices(right)])
    _, final = kernel_over_dvr(product_rows, min(n, m))
    assert final >= DIGITS-8


def operator(c, basis):
    columns = []
    for e in basis:
        image = (c.transpose(2, 3, 0, 1).reshape(9, 9)@e.T.flatten()).reshape(3, 3)
        columns.append([image[i, j]-int(i == j)*image[2, 2] for i, j in INDICES])
    return np.array(columns, dtype=object).T % MODULUS


def valuation(value):
    value = int(value) % MODULUS
    if not value:
        return DIGITS
    n = 0
    while value % P == 0:
        value //= P
        n += 1
    return n


def main():
    data = lift()
    basis = units(MODULUS)
    c = brackets(data[2], data[3]) % MODULUS
    check_inverse(c, basis)
    check_inverse(c.transpose(2, 3, 0, 1), basis)
    left, right = operator(c, basis), operator(c.transpose(2, 3, 0, 1), basis)
    char = sp.Matrix(right@left % MODULUS).charpoly().as_poly()
    poly = sp.Poly.from_list([int(v) % MODULUS for v in char.all_coeffs()], gens=char.gens)
    assert valuation(sp.Matrix(left).det()) == 8
    assert valuation(poly.discriminant()) == 62
    print('PASS: 62 pivots, both 18+18 trace ranks and 36 product ranks,'
          ' invertible crossed bracket and eight simple spectral lines')


if __name__ == '__main__':
    main()
