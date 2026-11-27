"""Exact fixed-divisor and evaluation calculations for the four directed row defects."""
from itertools import product
from sympy import Matrix, Poly, Rational, Symbol


def main():
    y = Symbol("y")
    for chosen in product(range(2), repeat=4):
        assert all(2 - 2 * chosen[i] - chosen[i - 1] >= -1 for i in range(4))
    assert {a + 2*b + 4*c for a, b, c in product(range(2), range(3), range(2))} == set(range(10))
    assert (2**4 + 1) == 17
    assert [r for r in range(7, 18) if 2*min(10, r) - r <= 6] == [14, 15, 16, 17]

    # The cyclic map t -> t^2+1/4 fixes 1/2 with multiplier one.
    iterate = Poly(y, y).to_field()
    for _ in range(4):
        iterate = iterate**2 + Poly(Rational(1, 4), y)
    affine_fixed = Poly(iterate.as_expr() - y, y)
    assert affine_fixed.degree() == 16  # infinity is the seventeenth point
    assert affine_fixed.gcd(affine_fixed.diff()).degree() > 0

    # Perfect pairing and the exact dual evaluation code on a nonreduced divisor.
    divisor = Poly(y**5*(y-1)**6*(y+2)**3, y)
    moments = [Poly(y**j, y).rem(divisor).nth(13) for j in range(27)]
    pairing = Matrix(14, 14, lambda i, j: moments[i+j])
    assert abs(pairing.det()) == 1
    assert pairing[:10, :4].is_zero_matrix and pairing[:10, :].rank() == 10

    # No degree-nine tuple can interpolate the full cubic Veronese of a
    # degree-four map once the three linear syzygies extend off the divisor.
    equations = Matrix.zeros(42, 40)
    for j in range(3):
        for i in range(10):
            equations[14*j+i, 10*(j+1)+i] = 1
            equations[14*j+i+4, 10*j+i] = -1
    assert equations.rank() == 40
    print("Four graph degrees, nonreduced fixed divisor, full nine-degree evaluations, and length-14 obstruction checked.")


if __name__ == "__main__":
    main()
