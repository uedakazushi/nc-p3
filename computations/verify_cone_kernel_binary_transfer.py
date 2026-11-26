"""Check the cone's cyclic commutator argument and its binary contraction bounds."""
from sympy import Matrix, Poly, Rational, eye, rem, symbols, zeros

t = symbols("t")
p = (t - 1)**4 * (t + 2)**2


def coordinates(f):
    poly = Poly(rem(f, p, t), t)
    return Matrix([poly.nth(j) for j in range(6)])


def multiplication(f):
    return Matrix.hstack(*(coordinates(t**j * f) for j in range(6)))


T = multiplication(t)
metric = Matrix(6, 6, lambda i, j: coordinates((1+t*t)*t**(i+j))[5])
assert metric.det() != 0 and metric*T == T.T*metric
u, v = coordinates(1+t), coordinates(2+t**3)
C = (u*v.T - v*u.T)*metric
assert C.rank() == 2 and C.T*metric == -metric*C
assert all((T**i*C).trace() == 0 for i in range(6))
assert (C*Matrix.hstack(*C.nullspace())).is_zero_matrix
# Solve the ordinary commutator system at a repeated-root cyclic operator.
def matrix_unit(i, j):
    a = zeros(6)
    a[i, j] = 1
    return a
units = [matrix_unit(i, j) for i in range(6) for j in range(6)]
operator = Matrix.hstack(*((e*T-T*e).reshape(36, 1) for e in units))
values, parameters = operator.gauss_jordan_solve(C.reshape(36, 1))
S = values.subs({s: 0 for s in parameters}).reshape(6, 6)
assert S*T - T*S == C and S.T*metric == metric*S

# A deficient binary model with dim K=2 and dim BK=4.
# Monomial vectors in S^d are indexed by their y exponent.
def mul_x(column):
    return column.col_join(zeros(1, 1))
def mul_y(column):
    return zeros(1, 1).col_join(column)
D = Matrix([[1,0,0,2], [0,1,0,0], [2,0,1,0], [0,2,0,0],
            [1,0,0,0], [0,1,0,1]])
K = Matrix.hstack(*(mul_x(D[:, i+1])-mul_y(D[:, i]) for i in range(3)))
BK = Matrix.hstack(*(m(K[:,i]) for m in (mul_x,mul_y) for i in range(3)))
assert K.rank() == 2 and BK.rank() == 4
G = Matrix.hstack(*BK.T.nullspace())
Jx, Jy = G[:7,:], G[1:,:]
assert Matrix.hstack(Jx,Jy).rank() <= 5
block = Jx.row_join(Jy).row_join(zeros(7,4)).col_join(
    zeros(7,4).row_join(Jx).row_join(Jy))
assert 12-block.rank() >= 2
print("PASS: repeated-root adjoints, commutator rank, binary apolar prolongation and path blocks")
