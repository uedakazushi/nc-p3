"""Check binary fiber balancing, splitting bounds, and the sextic genus budget."""
from itertools import product

from sympy import Matrix


def fiber_test(dual_numbers):
    # Bases (1,e) for k[e]/(e^2) or k[e]/(e^2-e), including the branch fiber.
    def multiply(i, j):
        if i == j == 1 and dual_numbers:
            return (0, 0)
        return (0, 1) if i or j else (1, 0)

    multiplication = Matrix.hstack(*(Matrix(multiply(i, j)) for i, j in product(range(2), repeat=2)))
    relations = multiplication.nullspace()
    words = list(product(range(2), repeat=3))
    left = Matrix([[v[2*i+j] if k == a else 0 for i, j, k in words]
                   for v in relations for a in range(2)])
    right = Matrix([[v[2*j+k] if i == a else 0 for i, j, k in words]
                    for v in relations for a in range(2)])
    assert len(relations) == 2 and left.rank() == right.rank() == 4
    assert left.col_join(right).rank() == 6


def main():
    for ramified in (False, True):
        fiber_test(ramified)
    for a in range(-1, 2):
        b = -2-a
        assert max(0, a+1) + max(0, b+1) <= 2
    assert 32-8-8+2 < 20
    assert [g for g in range(5) if 18-4*6+2*g >= 0] == [3, 4]
    print("PASS: both binary fibers, split intersection bounds, genus range and cubic exclusion")


if __name__ == "__main__":
    main()
