"""Verify the equality-case surface models and the finite-triple obstruction."""
from verify_cofactor_nineteen_clusters import blow_up_corners, fixed_points, intersection
from verify_cofactor_nineteen_surface import ruled_classes


def main():
    cone_models = ruled_classes(2, 2)
    surviving = [(a, b, base_length) for a, b, base_length in cone_models
                 if b - 2 * a <= 2]
    assert surviving == [(2, 5, 4)]
    smooth_models = [(a, b) for a, b, _ in ruled_classes(0, 1)
                     if a <= b and (a - 1) * (b - 1) == 2]
    assert smooth_models == [(2, 3)]
    general_length = intersection(blow_up_corners(fixed_points(2, 3, 11, 47)))
    assert general_length == 16 < 20 - 3
    print("PASS: genus-two quadric models, length four base scheme and at least seventeen isolated fibers in the finite case")


if __name__ == "__main__":
    main()
