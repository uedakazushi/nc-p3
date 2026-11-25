"""Symbolic PBW and change-of-lift computations for every nondiagonal plane quotient."""

import sympy as sp

ALT = []
for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
    form = sp.zeros(3)
    form[j, k], form[k, j] = 1, -1
    ALT.append(form)


def pbw(h, mat):
    value = (mat*(h*h).T-h.det()*h.inv()*mat.T)*h.inv().T
    return value, sp.Matrix([value[1, 2], value[2, 0], value[0, 1]])


def shifts(h):
    relations = [r*h.T for r in ALT]
    return sp.Matrix.hstack(*[sp.Matrix.hstack(*[(r+r.T)[:, i] for r in relations]).reshape(9, 1)
                              for i in range(3)])


def normal_space(h):
    v = sp.symbols('v:9')
    f, _ = pbw(h, sp.Matrix(3, 3, v))
    linear = (f+f.T).reshape(9, 1).jacobian(v)
    move = shifts(h)
    pivot = list(move.T.rref()[1])
    return linear.col_join(sp.eye(9)[pivot, :])


def act(h, g, mat, const):
    move = shifts(h)
    rows, cols = list(move.T.rref()[1]), list(move.rref()[1])
    raw = g*mat*g.T/g.det()
    e = sp.zeros(3, 1)
    for i, v in zip(cols, move[rows, cols].inv()*raw.reshape(9, 1)[rows, :]):
        e[i] = v
    relations = [r*h.T for r in ALT]
    normalized = raw-(move*e).reshape(3, 3)
    result = g*const/g.det()+raw.T*e-sp.Matrix([(e.T*r*e)[0] for r in relations])
    return normalized.applyfunc(sp.factor), result.applyfunc(sp.factor)


def verify():
    h = sp.eye(3)
    h[0, 1] = 1
    a, b, c, d = sp.symbols('a b c d')
    base = sp.Matrix([[0, -a/2, 0], [a, 0, -2*c], [b, c, d]])
    linear = normal_space(h)
    assert linear.rank() == 5 and linear*base.reshape(9, 1) == sp.zeros(linear.rows, 1)
    b0, b1, b2 = sp.symbols('b0 b1 b2')
    const = sp.Matrix([b0, b1, b2])
    _, coeff = pbw(h, base)
    equations = base*coeff+(h*h-h.inv())*const
    assert (equations-sp.Matrix([3*b1-a*(b-c)/2, 0, -(3*a*d+4*b*c+2*c*c)/2])).applyfunc(sp.expand).is_zero_matrix
    assert sp.expand(const.dot(coeff)+3*c*b0-(b-c)*b1+3*a*b2/2) == 0
    p, q, r, z, t = sp.symbols('p q r z t', nonzero=True)
    g = sp.Matrix([[p, q, r], [0, p, 0], [0, z, t]])
    transformed, residual = act(h, g, base, const.subs(b1, a*(b-c)/6))
    expected = [(a*p-2*c*r)/(p*t), (3*a*p*z+4*b*p*t-6*c*r*z+6*d*r*t)/(4*p*p*t), c/p, (d*t-c*z)/(p*p)]
    assert all(sp.factor(transformed[i]-val) == 0 for i, val in zip([3, 6, 7, 8], expected))
    assert residual.subs({a:0, b:0, c:0, d:1, r:0, t:p*p}).applyfunc(sp.factor) == sp.Matrix([b0/p**3, 0, b2/p**2])
    diff = residual.subs({a:0, b:1, c:0, d:0, p:1})-sp.Matrix([(b0+r*(b2-1))/t, 0, b2])
    assert diff.applyfunc(sp.factor).is_zero_matrix

    h[1, 2] = 1
    base = sp.Matrix([[0, 0, 0], [2*a/3, -a/3, 0], [a, 0, 0]])
    linear = normal_space(h)
    assert linear.rank() == 8 and (linear*base.reshape(9, 1)).is_zero_matrix
    _, coeff = pbw(h, base)
    rhs = sp.Matrix([b0, 0, a*a/9])
    assert (base*coeff+(h*h-h.inv())*rhs).applyfunc(sp.expand).is_zero_matrix and coeff.dot(rhs) == 0
    g = sp.Matrix([[p, q, r], [0, p, q], [0, 0, p]])
    m2, beta = act(h, g, base, rhs)
    assert (m2-base/p).applyfunc(sp.factor).is_zero_matrix
    assert sp.factor(beta[0].diff(r)+8*a*a/(9*p**3)) == 0

    h[1, 2], h[2, 2] = 0, t
    linear = normal_space(h)
    rows = list(linear.subs(t, 2).T.rref()[1])
    det = sp.factor(linear[rows, :].det())
    for expr in sp.factor_list(det.as_numer_denom()[0], t)[1]:
        assert expr[0] in (t, t-1, t+1)
    special = normal_space(h.subs(t, -1))
    assert special.nullspace() == [sp.eye(9)[:, 8]]


if __name__ == '__main__':
    verify()
    print('PASS: all nondiagonal Jordan types, symbolic PBW, residual scalings, exceptional split eigenvalues')
