"""Check the moving-plane tangent condition and the affine residual equations."""
from sympy import Matrix, diff, expand, symbols


def main():
    x, y, z = coordinates = symbols("x y z")
    p, q, s = symbols("p q s")
    upper = symbols("m:6")
    m = Matrix([[upper[0], upper[1], upper[2]],
                [upper[1], upper[3], upper[4]],
                [upper[2], upper[4], upper[5]]])
    b = Matrix(symbols("b:3"))
    v = Matrix(coordinates)
    phi = s*(x**3+y**3+z**3)/3+(p+q)*x*y*z-(v.T*m*v)[0]/2-(b.T*v)[0]
    residual = Matrix([s*x*x+(p+q)*y*z, s*y*y+(p+q)*z*x, s*z*z+(p+q)*x*y])-m*v-b
    assert [expand(diff(phi, t)-r) for t, r in zip(coordinates, residual)] == [0]*3
    # The two normal directions in the central chart are the two point columns.
    u = Matrix(symbols("u:3"))
    mixed = u.row_join(-v)
    cross = u.cross(v)
    assert [mixed.extract(rows, [0, 1]).det() for rows in [(1, 2), (0, 2), (0, 1)]] == [
        -cross[0], cross[1], -cross[2]]
    print("PASS: transverse minors and all three critical-scheme equations")


if __name__ == "__main__":
    main()
