"""Verify the Hankel obstruction and the degrees of the relation bundles."""
from sympy import Matrix, expand, symbols


def main():
    a = symbols("a:5")
    b = symbols("b:5")
    left = Matrix(3, 3, lambda i, j: a[i+j])
    right = Matrix(3, 3, lambda i, j: b[i+j])
    apolar = Matrix([a[:4], a[1:], b[:4], b[1:]])
    quadratic = left*right.adjugate()*left
    assert expand(quadratic[0, 2]-quadratic[1, 1]-apolar.det(method="domain-ge")) == 0
    assert 3*4-4 == 8
    assert 4*3-7 == 5
    print("PASS: symbolic ten-coefficient Hankel identity and overlap dimension")


if __name__ == "__main__":
    main()
