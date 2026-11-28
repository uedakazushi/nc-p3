"""Exact weighted path-presentation bounds for generic simple extensions."""

from sympy import Matrix, expand, symbols


def main():
    d = symbols("d:4")
    euler = sum(d[i] * (2*d[i] - 4*d[i-1] - 4*d[(i+1)%4] + 6*d[(i+2)%4])
                for i in range(4))
    signature = 4*(d[0]+d[2]-d[1]-d[3])**2 - 2*((d[0]-d[2])**2+(d[1]-d[3])**2)
    assert expand(euler-signature) == 0

    # Equitable row sums for the two vertex types, with any number of middle vertices.
    m, c, epsilon, p, a_sum, b_sum, b, r, k = symbols("m c epsilon p a_sum b_sum b r k")
    n, t = m**2, c/m
    middle_sum = m*(c*b_sum-r*epsilon)
    tip = (1 - t*m*(c*p-epsilon*a_sum) + t**2*(n*(p-1)+1)
           + k*t**3*middle_sum + k*t**4*n)
    center = m*(c*b-epsilon) - t*n*b + k*t**2*middle_sum + k*t**3*n
    assert expand(tip).coeff(m, 0) == 1-c**2+c*epsilon*a_sum
    assert expand(center).coeff(m, 1) == -epsilon
    assert expand(center).coeff(m, 0) == 0
    assert expand(m**2*(tip-(1-c**2+c*epsilon*a_sum))).has(m) is False
    assert expand(m*(center+m*epsilon)).has(m) is False

    # Kernel of a rank-factorized tensor as relations on two-step words.
    left, right = Matrix([[1, 0], [2, 1], [0, 1]]), Matrix([[1, 2, 1], [0, -1, 3]])
    tau = left*right
    assert tau.rank() == 2
    assert Matrix(1, 9, list(tau)).rank() == 1
    print("Euler signature and both weighted growth bounds verified symbolically.")


if __name__ == "__main__":
    main()
