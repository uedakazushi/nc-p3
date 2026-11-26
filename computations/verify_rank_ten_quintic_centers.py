"""Determinant identity for a square factor at the unique quintic center."""
import sympy as sp


def main():
    s, t, e, f, b1, b2, h = sp.symbols("s t e f b1 b2 h")
    linear = sp.Matrix([e, f, 0])
    fixed = sp.Matrix([0, 0, 1])
    last = sp.Matrix([b1, b2, s*h])
    matrix = sp.Matrix.hstack(s*linear, t*linear, s*s*fixed, last)
    residual = e*b2-f*b1
    assert sp.expand(matrix[:, (0, 2, 3)].det()+s**3*residual) == 0
    assert sp.expand(matrix[:, (1, 2, 3)].det()+t*s**2*residual) == 0
    for minor in (matrix.subs(s, 0)[:, (1, 3)],):
        assert sp.expand(minor[:2, :].det()-t*residual) == 0
        assert minor[2, :] == sp.zeros(1, 2)
    print("PASS: concentrated quintic and rank drop after removing the square factor")


if __name__ == "__main__":
    main()
