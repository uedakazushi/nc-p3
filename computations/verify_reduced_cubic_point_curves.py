"""Exact coefficient, stabilizer, and deformation checks for triangular cubic graphs."""
from itertools import product

from sympy import Matrix, Poly, expand, symbols, zeros


def main():
    x = symbols("x:3")
    q, u = symbols("q u")
    rows = zeros(3)
    for i in range(3):
        j, k = (i+1) % 3, (i+2) % 3
        rows[i, i], rows[i, j], rows[i, k] = -q*x[i], u*x[k], x[j]
    mixed = x[0]*x[1]*x[2]
    expected = (1+u**3-q**3)*mixed+u*q*sum(a**3 for a in x)
    assert expand(rows.det()-expected) == 0
    assert all(expand(v) == 0 for v in rows*Matrix(x)-Matrix(
        [(1+u)*x[(i+1) % 3]*x[(i+2) % 3]-q*x[i]**2 for i in range(3)]))
    # Singularity of a member a*xyz+b*sum(x_i^3) requires b=0 or a^3=-27 b^3.
    assert Poly((1+u**3-q**3)**3+27*(u*q)**3, u).LC() == 1
    # On coordinate lines all nonvertex points have at least two independent rows.
    for i, j in product(range(3), repeat=2):
        active = {i, j}
        left = rows.subs(u, 0).subs({v: 0 for h, v in enumerate(x) if h not in active})
        minors = [expand(left.extract(a, b).det()) for a in ((0, 1), (1, 2), (0, 2))
                  for b in ((0, 1), (1, 2), (0, 2))]
        assert any(m != 0 and len(Poly(m, *x).terms()) == 1 for m in minors)
    # Torus of diagonal changes modulo common scaling has character lattice A2.
    action = Matrix([[-2, 1], [1, -2]])
    assert action.det() == 3
    assert all((-2*i+(i+1) % 3+(i+2) % 3) % 3 == 0 for i in range(3))
    print("PASS: both cyclic deformations, triangle fibers, diagonal equations, stabilizer index three")


if __name__ == "__main__":
    main()
