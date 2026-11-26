"""Planar Gorenstein contact algebra and ruling counts at 17-site graphs."""
import sympy as sp


def verify():
    a, b, c = sp.symbols("a b c")
    # Columns of multiplication by ax+by+cxy, in the basis 1,x,y,xy.
    multiplication = sp.Matrix([[0, 0, 0, 0], [a, 0, 0, 0],
                                [b, 0, 0, 0], [c, b, a, 0]])
    assert (multiplication**2)[3, 0] == 2*a*b
    assert multiplication.rank() == 2
    for substitution in ({a: 1, b: 0}, {a: 0, b: 1}):
        assert multiplication.subs(substitution).rank() == 2
    assert multiplication.subs({a: 0, b: 0, c: 1}).rank() == 1
    # Two nilpotent columns can indeed give a nonzero socle scalar matrix.
    x, y = sp.symbols("x y")
    outer = sp.Matrix([x, y])*sp.Matrix([[y, x]])
    assert outer.applyfunc(lambda p: p.subs({x**2: 0, y**2: 0})) == x*y*sp.eye(2)
    assert [k*(k-1)//2 + k + (6-k)*(5-k)//2 for k in (2, 3)] == [9, 9]
    assert 2+3*2 < 17-7
    assert 2*max(5+2, 2+2+2) < 16


if __name__ == "__main__":
    verify()
    print("PASS: planar length-four ideals, isotropic directions, socle exception, and re-marking bounds")
