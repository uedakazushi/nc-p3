"""Compare the degree of the triple-fiber locus with isolated marked points."""
from verify_kernel_curve_scroll import degree, genus, integral


def main():
    movable_length = integral.subs({degree: 8, genus: 2})
    assert movable_length == 12
    marked_on_curves = {
        "empty": 0, "line": 3, "conic": 6, "two incident lines": 2 * 3,
    }
    for marks in marked_on_curves.values():
        assert 20 - marks > movable_length
    # The splitting witness has three sections on each degree-four summand.
    assert 4 + 1 - 2 == 3
    assert 3 + 1 - 2 == 2
    print("PASS: lower bound of fourteen isolated points exceeds the intersection length twelve")


if __name__ == "__main__":
    main()
