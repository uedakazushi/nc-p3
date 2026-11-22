"""Degree-six ideal ranks from saturated degree-five central products on the trace disk."""

from itertools import product

from verify_klein_center_lift import add_multiple, echelon
from verify_rank_fifteen_cubic_zero_disk import pivot_rows
from verify_rank_fifteen_quadratic_dual import sixth_insertion
from verify_rank_fifteen_variable_zero_disk import Jet, center_data
from verify_rank_fourteen_all_ordinary import saturated_projection
from verify_rank_fourteen_geometric_reduction import Scalar


def verify():
    windows, basis, image, centers = center_data(9)
    precision = Jet.precision
    product_spaces = {}
    for start, h in product(range(4), repeat=2):
        rows = []
        for g, center in centers[start]:
            row = {}
            for word, coefficient in center.items():
                add_multiple(row, image(start, word+(g^h,)), coefficient)
            rows.append(row)
        width = len(basis[start, 5, h])
        pivots, loss = pivot_rows(rows, width, 10, Jet)
        assert width == 14 and loss == 4
        product_spaces[start, h] = list(pivots.values())

    for start, h in product(range(4), repeat=2):
        Jet.precision = precision
        rows, width = sixth_insertion(windows, basis, image, start, h)
        _, projection, loss = saturated_projection(rows, width, (20, 20, 22, 22)[h], Jet)
        assert loss == 2
        # The saturated degree-five bases retain only the constant coefficient.
        Jet.precision = 1
        residues = []
        offset = 0
        for last in range(4):
            for row in product_spaces[start, h^last]:
                projected = {}
                for j, coefficient in row.items():
                    add_multiple(projected, projection[offset+j], coefficient)
                residues.append({j: Scalar(c.cs[0].specialize())
                                 for j, c in projected.items() if c.unit()})
            offset += len(basis[start, 5, h^last])
        assert offset == width and len(residues) == 16
        assert len(echelon(residues)) == (8, 8, 10, 10)[h]
    print("All degree-six residue product ranks are (8,8,10,10); quotient dimension is 44 to 48.")


if __name__ == "__main__":
    verify()
