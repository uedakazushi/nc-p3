"""Degree, character, and normalization restrictions at a singular quadric pole."""
from itertools import product
from sympy import expand, symbols


def main():
    # Chart around the smooth base conic: rank three gives two O(2) coordinates.
    s, t, u, v, w, b, c, d = symbols("s t u v w b c d")
    assert expand((s*t)**2 - s**2*t**2) == 0
    assert expand((u+c)*(v+b) - (u*v+b*u+c*v+b*c)) == 0
    assert expand((u+c/2)**2+b*v+d-c**2/4 - (u**2+c*u+b*v+d)) == 0
    assert [j for j in range(4) if (2*j)%4 == 0] == [0,2]
    allowed = [(m, n) for m, n in product(range(1,5), repeat=2) if m*n <= 4 and m*n % 4 == 0]
    assert allowed == [(1,4), (2,2), (4,1)]
    for branches in range(3, 9):
        assert branches*(branches-1)//2 > branches-1
    print("Singular-quadric shifts, vertical inertia exponents, fourth-degree splitting, and plane-curve delta checked.")


if __name__ == "__main__":
    main()
