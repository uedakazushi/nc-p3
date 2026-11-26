"""Verify the Euler characteristic identities and Hilbert function bounds."""
from math import comb

from sympy import expand, solve, symbols


def hilbert_genus_bound(degree, first_differences):
    """Lower bounds for the general linear quotient, constant thereafter."""
    return sum(degree - h for h in first_differences[1:])


def main():
    d, g, z, n = symbols("d g z n")
    # Five (1,1) equations in P^3 x P^3. Evaluate the polynomial, including
    # negative twists, on the diagonal and both coordinate axes.
    def chi(t):
        return (t + 1) * (t + 2) * (t + 3) / 6

    a, b = symbols("a b")
    chi_w = expand(sum((-1)**i * comb(5, i) * chi(a-i) * chi(b-i)
                       for i in range(6)))
    assert chi_w == 10*a + 10*b - 10
    chi_residual_ideal = z - 1 - 2*(10-d)
    assert solve(chi_w.subs({a: 0, b: 0}) - chi_residual_ideal - (1-g), z) == [g + 10 - 2*d]
    defect = expand(2*(g + 10 - 2*d) - 2)
    assert defect == 18 - 4*d + 2*g
    assert hilbert_genus_bound(4, [1, 3, 4]) == 1
    assert hilbert_genus_bound(5, [1, 3, 3, 4, 5]) == 5  # coarse bound
    assert hilbert_genus_bound(5, [1, 3, 6, 4, 5]) == 2
    plane_genus = lambda e: (e-1)*(e-2)//2
    for e in (1, 2):
        total = plane_genus(e) + plane_genus(5-e) + e-1
        assert total == (3 if e == 1 else 2)
    assert [(i, int(defect.subs({g: i, d: 5}))) for i in (1, 2)] == [(1, 0), (2, 2)]
    assert defect.subs({g: 3, d: 6}) == defect.subs({g: 5, d: 7}) == 0
    # Euler characteristic of the conormal presentation. Its kernel is the
    # dual of the fixed ideal, with twists -H-L2-L(-1) of total degree -4d.
    f, j = symbols("f j")
    degree_source, degree_target = -12*d, -8*d-(2*g-2)
    assert solve(degree_source - degree_target + j - (-4*d + f), f) == [2*g+j-2]
    print("PASS: linkage and conormal polynomials, genus bounds, and exceptional lengths")


if __name__ == "__main__":
    main()
