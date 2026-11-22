"""Independence of all sixteen central-quartic times arrow products on the trace disk."""

from itertools import product

from verify_klein_center_lift import add_multiple
from verify_rank_fifteen_variable_zero_disk import Jet, center_data
from verify_rank_fourteen_all_ordinary import saturated_projection


def verify():
    _, basis, image, centers = center_data(9)
    assert Jet.precision == 5
    for start, h in product(range(4), repeat=2):
        rows = []
        for g, center in centers[start]:
            arrow = g ^ h
            row = {}
            for word, coefficient in center.items():
                add_multiple(row, image(start, word+(arrow,)), coefficient)
            rows.append(row)
        width = len(basis[start, 5, h])
        assert (width, len(rows)) == (14, 4)
        _, _, loss = saturated_projection(rows, width, 10, Jet)
        assert loss == 4
    print("All 16 degree-five product blocks have rank four and quotient dimension ten.")


if __name__ == "__main__":
    verify()
