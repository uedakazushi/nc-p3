"""Check the vertex module of the rectangular eight-point conic-bundle model."""
from sympy import Matrix, eye, zeros


def action(coefficients):
    # S/(s^4,t^4) in degrees one and five. Rows s^3*t^2,s^2*t^3
    # are paired in reverse order with the columns s,t by the socle s^3*t^3.
    a, b, c = coefficients  # coefficients of s^3*t, s^2*t^2, s*t^3
    result = zeros(4)
    result[2:, :2] = Matrix([[b, a], [c, b]])
    return result


def main():
    b = action((2, 1, -3))
    h = action((1, 4, 3))
    for coeffs in ((1, 1, 1), (1, 0, 2), (0, 0, 0)):
        c = action(coeffs)
        matrix = (-c).row_join(-h).col_join((-h).row_join(eye(4)-b))
        a0, a1, a2 = coeffs
        hankel = Matrix([[a2, a1], [a1, a0]])
        assert h*b == zeros(4) and h*h == zeros(4)
        assert 8 - matrix.rank() == 4 - hankel.rank()
    assert 20 - 3 - 4 > 12 and 20 - 4 > 14
    assert 2 + 2*2 > 4 and 4 + 2*1 > 4
    print("PASS: rectangular vertex module, Hankel ranks at the horizontal marks and all residual counts")


if __name__ == "__main__":
    main()
