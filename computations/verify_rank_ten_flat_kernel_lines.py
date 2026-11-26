"""Check the local minors detecting the limiting successor line."""
import itertools

import sympy as sp


def main():
    t, c1, c2, c3, c4 = sp.symbols("t c1 c2 c3 c4")
    matrix = sp.Matrix([[1, 0, 0, 0], [0, t, 0, 0], [c1, c2, c3, c4]])
    minors = [matrix[:, idx].det() for idx in itertools.combinations(range(4), 3)]
    assert minors == [t*c3, t*c4, 0, 0]
    a, b, x, y, s, v = sp.symbols("a b x y s v")
    u, v0, v1 = (sp.eye(3).col(i) for i in range(3))
    pencil = sp.Matrix.hstack(s*s*u, s*v*u, s*(s*v0+v*v1), v*(s*v0+v*v1))
    point = sp.Matrix([-b*x, a*x, -b*y, a*y])
    value = pencil * point
    expected = (a*v-b*s)*(s*x*u+y*(s*v0+v*v1))
    assert (value-expected).applyfunc(sp.expand) == sp.zeros(3, 1)
    print("PASS: local deficit minors and planar evaluation on every limiting line")


if __name__ == "__main__":
    main()
