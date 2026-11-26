"""Bounds and local identities for a graph with one dual-number factor."""
import sympy as sp


def main():
    assert 19-7 == 12 and 20-7 == 13
    assert max(2+3+6, 6+4+1, 7+4) < 12
    assert 6+3*2 < 13 and (18+1)//2 > 7
    x, a, b = sp.symbols("x a b")
    pencil = sp.Matrix([[x, 0, a], [0, x, b]])
    assert [pencil[:, cols].det() for cols in ((0, 1), (0, 2), (1, 2))] == [x*x, x*b, -x*a]
    s, t, h = sp.symbols("s t h")
    col = sp.Matrix([s, t, h])
    rank_one = col*col.T
    tangent = rank_one.diff(h).subs(h, 0)
    assert tangent[:2, 2] == sp.Matrix([s, t]) and tangent[2, 2] == 0
    print("PASS: nineteen-site counts, projected ruling minors and transverse conic tangent")


if __name__ == "__main__":
    main()
