"""Exact determinantal identities for the pencils at a ternary branch center."""
from sympy import Matrix, eye, symbols


def main():
    s, t, a, b, c = symbols("s t a b c")
    pencil = Matrix([[s, t, 0], [0, 0, s], [0, 0, t]])
    coeff = [pencil.diff(v) for v in (s, t)]
    last = Matrix.vstack(*[(d*Matrix([a, b, c])).T for d in coeff])
    minors = [last[:, [i, j]].det() for i, j in ((0, 1), (0, 2), (1, 2))]
    assert all(f.subs(c, 0) == 0 for f in minors)
    assert c**2 in minors
    assert last.subs(c, 0) == Matrix([[a, 0, 0], [b, 0, 0]])
    jordan = Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    for block in (jordan, jordan + eye(3), Matrix.diag(1, 2, 3)):
        assert (s*eye(3) + t*block).det().coeff(s, 3) == 1
    print("PASS: singular ternary pencil, fixed pure factor, regular determinant examples")


if __name__ == "__main__":
    main()
