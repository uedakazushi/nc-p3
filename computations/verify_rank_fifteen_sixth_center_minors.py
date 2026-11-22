"""Two H divisions and a unit minor in every unresolved sixth-degree character."""

from itertools import combinations, product

from flint import fmpz_poly

from verify_klein_center_lift import add_multiple
from verify_rank_fifteen_cubic_zero_disk import pivot_rows
from verify_rank_fifteen_generic_levels import CommonDenominatorJet
from verify_rank_fifteen_quadratic_dual import sixth_insertion
from verify_rank_fifteen_variable_zero_disk import center_data, trace_relations
from verify_rank_fourteen_all_ordinary import saturated_projection

ONE = fmpz_poly([1])
MODULUS = fmpz_poly([2, 1, 1])


def residue_nonzero(poly):
    remainder = poly % MODULUS
    return bool(int(remainder[0]) % 3 or int(remainder[1]) % 3)


class TraceJet(CommonDenominatorJet):
    """H-jets with a shared Z[z] denominator, localized at (3,z²+z+2)."""

    def unit(self):
        assert residue_nonzero(self.den)
        return bool(self.coeffs) and residue_nonzero(self.coeffs[0])


def relations():
    one = TraceJet.constant(1)
    half = TraceJet([ONE], fmpz_poly([2]))
    z = TraceJet([fmpz_poly([0, 1])])
    h = TraceJet([fmpz_poly([]), ONE])
    return trace_relations(one, half, z, h, TraceJet.precision)


def remove_unit_pivots(rows):
    """Eliminate over the local ring, retaining original row and column labels."""
    residual = [(i, row.copy()) for i, row in enumerate(rows)]
    count = 0
    while True:
        pivot = next(((i, k)
                      for i in sorted(range(len(residual)), key=lambda j: len(residual[j][1]))
                      for k in sorted(residual[i][1]) if residual[i][1][k].unit()), None)
        if pivot is None:
            return count, residual
        i, k = pivot
        _, row = residual.pop(i)
        inverse = row.pop(k).inverse()
        normalized = {j: c*inverse for j, c in row.items()}
        for _, other in residual:
            if k in other:
                add_multiple(other, normalized, -other.pop(k))
        count += 1


def unit_minor(residual):
    zero = TraceJet.constant(0)
    for (i, a), (j, b) in combinations(residual, 2):
        for k, l in combinations(sorted(a.keys() | b.keys()), 2):
            det = a.get(k, zero)*b.get(l, zero) - a.get(l, zero)*b.get(k, zero)
            if det.unit():
                return i, j, k, l
    raise AssertionError("Missing rank-two unit minor")


def verify():
    windows, basis, image, centers = center_data(11, TraceJet, relations)
    precision = TraceJet.precision
    assert precision == 7
    for start in range(4):
        TraceJet.precision = precision
        spaces = {}
        for h in range(4):
            rows = []
            for g, center in centers[start]:
                row = {}
                for word, coefficient in center.items():
                    add_multiple(row, image(start, word+(g^h,)), coefficient)
                rows.append(row)
            pivots, loss = pivot_rows(rows, len(basis[start, 5, h]), 10, TraceJet)
            assert loss == 4
            spaces[h] = list(pivots.values())

        for h in (2, 3):
            TraceJet.precision = precision
            rows, width = sixth_insertion(windows, basis, image, start, h)
            _, projection, loss = saturated_projection(rows, width, 22, TraceJet)
            assert width == 56 and loss == 2
            TraceJet.precision = precision-4
            appended = []
            for last in range(4):
                for row in spaces[h^last]:
                    result = {}
                    for j, coefficient in row.items():
                        add_multiple(result, projection[14*last+j], coefficient)
                    appended.append(result)
            count, residual = remove_unit_pivots(appended)
            assert count == 10 and len(residual) == 6
            # Each call checks exact vanishing over Q(z) before dropping H.
            for _ in range(2):
                residual = [(i, {j: v.divide_by_parameter() for j, v in row.items()})
                            for i, row in residual]
                TraceJet.precision -= 1
            assert TraceJet.precision == 1
            witness = unit_minor(residual)
            print("Start", start, "character", h, "unit minor", witness, flush=True)


if __name__ == "__main__":
    verify()
