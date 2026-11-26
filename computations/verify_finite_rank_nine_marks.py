"""Exact checks for the linear jets and symmetric slices at a planar mark."""
from sympy import Matrix, Rational, eye, symbols, zeros


def main():
    x, y, s, n = symbols("x y s n")
    v = Matrix([x, y])
    coefficients = symbols("c:8")
    block = Matrix(2, 2, coefficients[:4])*x + Matrix(2, 2, coefficients[4:])*y
    # Independent row/column coordinates normalize both tangent isomorphisms.
    identities = list(block*v) + list(v.T*block)
    rows = [expr.expand().coeff(x, a).coeff(y, 2-a)
            for expr in identities for a in range(3)]
    assert Matrix(rows).jacobian(coefficients).rank() == 8

    for k in (Matrix([[1, 1], [-1, 0]]), Matrix([[1, 2], [0, 1]]), eye(2)):
        assert k.det() != 0
        a = Matrix([[s]]).row_join(-v.T*k).col_join((-k*v).row_join(n*k))
        kernel = Matrix([n, x, y])
        quadric = s*n - (v.T*k*v)[0]
        expected = Matrix([quadric, 0, 0])
        assert (a*kernel-expected).applyfunc(lambda e: e.expand()) == zeros(3, 1)
        assert (a.T*kernel-expected).applyfunc(lambda e: e.expand()) == zeros(3, 1)
        assert (a.det()-k.det()*n*quadric).expand() == 0
        symmetric = (a+a.T)*Rational(1, 2)
        image = Matrix(list(symmetric)).jacobian((s, x, y, n))
        assert image.rank() == 2 + (k+k.T).rank() >= 3
    print("PASS: two-sided first jets, common quadric kernels and symmetric-slice ranks")


if __name__ == "__main__":
    main()
