"""Check the binary sextic elimination and the nonreduced fixed intersection."""
from sympy import Matrix, factor, resultant, symbols


def main():
    d, b, u = symbols("d b u")
    reciprocals = [1/(1+j*d) for j in range(4)]
    numerators = [(reciprocals[j+2]-2*b*reciprocals[j+1]+b*b*reciprocals[j])
                  .together().as_numer_denom()[0] for j in (0, 1)]
    assert factor(resultant(*numerators, b)) == 8*d**6*(d+1)*(2*d+1)

    # Restriction of the two brackets to a local intersection of multiplicity
    # m0+m1. Multiplication by the fixed-line equation has kernel (u^m1).
    for m0 in range(3):
        for m1 in range(3):
            length = m0 + m1
            if not length:
                continue
            product = Matrix(length, length, lambda row, col: int(row == col+m0))
            assert product.rank() == m1
            assert all(product[:, j].is_zero_matrix for j in range(m1, length))
    # A Mobius transformation fixing a double point with first derivative -1
    # is an involution for all higher coefficients.
    z, c = symbols("z c")
    g = -z/(1+c*z)
    assert factor(g.subs(z, g)-z) == 0
    x = symbols("x:4")
    # Evaluation of (1,1) sections on a common length-four reduced frame.
    # Three normalized nodes are 0, 1, infinity in both copies.
    evaluation = Matrix([[1, 0, 0, 0], [1, 1, 1, 1], [0, 0, 0, 1],
                         [1, x[0], x[1], x[0]*x[1]]])
    assert factor(evaluation.det()) in (x[0]-x[1], x[1]-x[0])
    print("PASS: reciprocal tangent cubics, nonreduced annihilators, involutions and cross evaluation")


if __name__ == "__main__":
    main()
