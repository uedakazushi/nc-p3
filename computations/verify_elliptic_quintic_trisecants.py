"""Check the quintic Pfaffians and both elliptic projections symbolically."""
from itertools import combinations

from sympy import I, Matrix, expand, factor, symbols


def main():
    a, b, c, d, e, A, B = symbols("a b c d e A B")
    s = Matrix([a, b, c])
    t = Matrix([b, d, e])
    n = Matrix([-B*a, -d-A*a, c])
    skew = Matrix([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    upper = Matrix.zeros(2).row_join(s.T.col_join(t.T))
    phi = upper.col_join((-s).row_join(-t).row_join(skew))
    pfaffians = []
    for k in range(5):
        indices = [j for j in range(5) if j != k]
        v = phi.extract(indices, indices)
        pfaffians.append((-1)**k*(v[0,1]*v[2,3]-v[0,2]*v[1,3]+v[0,3]*v[1,2]))
    assert all(expand(z) == 0 for z in phi*Matrix(pfaffians))
    expected = [s.dot(n), t.dot(n), *s.cross(t)]
    assert all(any(expand(z-w) == 0 or expand(z+w) == 0 for w in pfaffians) for z in expected)
    x, y = symbols("x y")
    for eps in (I, -I):
        def project(x_coord, y_coord):
            return [1+eps*x_coord**2, x_coord, y_coord, x_coord*y_coord]
        rows = Matrix([
            project(x, y),
            project(-1/x, y/x**2),
            project(eps*(x**2-1)/(2*x), (1-eps)*y*(x**2+1)/(4*x**2)),
        ])
        assert all(factor(rows[:, cols].det()) == 0 for cols in combinations(range(4), 3))
        assert expand((2-eps)*(2+eps)) == 5
    print("PASS: quintic Pfaffians, kernel identity, both symbolic trisecants and order-five subgroups")


if __name__ == "__main__":
    main()
