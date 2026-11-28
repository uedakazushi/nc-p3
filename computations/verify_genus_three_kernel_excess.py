"""Check genus-three generating-section moduli dimensions and the residual length."""
from verify_kernel_curve_scroll import degree, genus, integral


def main():
    g = 3
    hom_dimension = 32 + 4 * (1 - g)
    assert hom_dimension == 24
    special = [3 + 1 + 4*5 - 1, 2 + 4*5 - 1, 4*6 - 4]
    strictly_semistable = [6 + 1 + 16 - 1, 3 + 2 + 16 - 1]
    assert max(special + strictly_semistable) < hom_dimension
    ext_dimension = -(3 - 5 + 1 - g)
    marked = 3 + 3 + ext_dimension - 1 + 16 - 1
    assert marked == hom_dimension and marked - 2 < hom_dimension
    assert integral.subs({degree: 8, genus: g}) == 8 < 20 - 6
    print("PASS: all special loci have dimension below 24, marked maximal sublines have dimension at most 24, triple length eight")


if __name__ == "__main__":
    main()
