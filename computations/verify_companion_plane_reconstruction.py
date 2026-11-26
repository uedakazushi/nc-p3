"""Exact kernel and cubic-pencil calculations, including a length-four local base."""
from itertools import combinations
import sympy as sp
from verify_prescribed_discriminant_seeds import bracket, cross


def check(metric, quadrics):
    v = sp.Matrix(sp.symbols('x y z'))
    basis = []
    for i in range(3):
        for j in range(i, 3):
            h = sp.zeros(3)
            h[i, j] = h[j, i] = 1
            basis.append(metric.inv()*h)
    forms = [(v.T*n*v)[0] for n in quadrics]
    normals = [metric.inv()*n for n in quadrics]
    assert all(n.trace() == 0 for n in normals)
    constraints = sp.Matrix([[h.trace() for h in basis]] +
                            [[(n*h).trace() for h in basis] for n in normals])
    columns = constraints.nullspace()
    for j in range(6):
        candidate = sp.eye(6)[:, j]
        if sp.Matrix.hstack(*columns, candidate).rank() > len(columns):
            # Extend within trace zero, projecting scalar matrices away afterwards.
            columns.append(candidate)
    block = [sum((c*h for c, h in zip(col, basis)), sp.zeros(3)) for col in columns]
    block = [b-sp.eye(3)*b.trace()/3 for b in block]
    pivots = sp.Matrix.hstack(*[b.reshape(9, 1) for b in block]).rref()[1]
    block = [block[j] for j in pivots]
    assert len(block) == 5
    frame = sp.Matrix.hstack(*[b.reshape(9, 1) for b in block])
    # [v]_x as an endomorphism for this metric; the common volume scalar is immaterial.
    skew = metric.inv()*cross(v)
    a = sp.Matrix([[(skew*bracket(g, h)).trace() for h in block] for g in block[:3]])
    x = sp.Matrix([a[1, 2], a[2, 0], a[0, 1], 0, 0])
    rad = v*v.T*metric - (v.T*metric*v)[0]*sp.eye(3)/3
    rad = frame.gauss_jordan_solve(rad.reshape(9, 1))[0]
    assert x[:3, :].jacobian(v).det() != 0
    assert (a*x).applyfunc(sp.expand).is_zero_matrix and (a*rad).applyfunc(sp.expand).is_zero_matrix

    minors, wedge = [], []
    for i, j in combinations(range(5), 2):
        rest = [k for k in range(5) if k not in (i, j)]
        minors.append((-1)**(i+j-1)*a[:, rest].det())
        wedge.append(x[i]*rad[j]-x[j]*rad[i])
    pivot = next(i for i in range(10) if sp.expand(wedge[i]) != 0)
    ratio = sp.cancel(minors[pivot]/wedge[pivot])
    assert ratio.is_Rational and ratio != 0
    assert all(sp.expand(l-ratio*r) == 0 for l, r in zip(minors, wedge))

    monomials = [v[0]**i*v[1]**j*v[2]**(3-i-j) for i in range(4) for j in range(4-i)]
    def coeff(f):
        p = sp.Poly(f, *v)
        return sp.Matrix([p.coeff_monomial(m) for m in monomials])
    coef = sp.Matrix.hstack(*[coeff(q*z) for q in forms for z in v])
    section = x[:3, :].jacobian(v).inv()
    lifts = section*rad[:3, :]
    val = sp.zeros(6, 3)
    for i, (j, k) in enumerate(((1, 2), (2, 0), (0, 1))):
        val[:, i] = coef.gauss_jordan_solve(coeff(v[j]*lifts[k]-v[k]*lifts[j]))[0]
    first, second = [(v.T*val[j:j+3, :]*v)[0].expand() for j in (0, 3)]
    scalar = sp.cancel(second/forms[0])
    assert scalar.is_Rational and scalar != 0 and sp.expand(first+scalar*forms[1]) == 0
    return ratio, scalar


if __name__ == '__main__':
    print("split", check(sp.eye(3), [sp.diag(1, -1, 0), sp.diag(0, 1, -1)]))
    metric = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    # x^2, y^2+xz.  At z=1 the base algebra is Q[y]/(y^4).
    print("curvilinear", check(metric, [sp.diag(1, 0, 0),
          sp.Matrix([[0, 0, sp.Rational(1, 2)], [0, 1, 0], [sp.Rational(1, 2), 0, 0]])]))
    print("PASS: polynomial kernel bivector and nonzero reconstruction scalar in both examples")
