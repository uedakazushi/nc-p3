"""The x=-1 residue direction at the pole of the integral-slope trace parameter."""

from itertools import product

from flint import fmpz_mod_poly_ctx, fmpz_poly

from verify_klein_center_lift import central_equations, field_paths
from verify_rank_fifteen_cubic_zero_disk import UNCUT, pivot_rows, tensor_relations
from verify_rank_fifteen_overlap_annulus import TraceJet, ONE, ZERO, poly
from verify_rank_fifteen_parameter_orbits import pattern
from verify_rank_fifteen_sixth_center_minors import TraceJet as RationalTraceJet
from verify_rank_fifteen_two_sided_cut import initial_products, overlap
from verify_rank_fourteen_all_ordinary import saturated_projection

MODULUS = poly([1, 1, 0, 2, 1])


class PoleJet(TraceJet):
    def unit(self):
        assert self.den % MODULUS
        return bool(self.coeffs) and bool(self.coeffs[0] % MODULUS)


class IntegralPoleJet(RationalTraceJet):
    modulus = fmpz_poly([1, 1, 0, 2, 1])

    def unit(self):
        def nonzero(poly):
            return any(int(c) % 3 for c in poly % self.modulus)
        assert nonzero(self.den)
        return bool(self.coeffs) and nonzero(self.coeffs[0])


def relations(jet, half, r, s):
    one = jet.constant(1)
    den = (one+s*s).inverse()
    a, b = (one-s*s)*den, 2*s*den
    t = a*b

    def root(square, initial):
        value = initial
        for _ in range(jet.precision.bit_length()):
            value = (value+square*value.inverse())*half
        assert not (value*value-square)
        return value

    k = root(one+r*r, one)
    p0_squared = 2*(k+one)
    p, u = (root(p0_squared+sign*2*t*r*r, jet.constant(2)) for sign in (1, -1))
    q, v = (r*root(4*p0_squared.inverse()+sign*2*t, initial)
            for sign, initial in ((1, a+b), (-1, b-a)))
    gamma = pattern(one, (k-1)*half, r*a*half, r*b*half, (u-p)*half*half,
                    (u+p)*half*half, (v-q)*half*half, (v+q)*half*half)
    return tensor_relations(gamma, one*0)


def initial_tables(jet, precision, half, r, s):
    jet.precision = precision
    windows = relations(jet, half, r, s)
    normalization = []
    order = iter(product(range(1, 6), range(4), range(4)))
    losses = []

    def reduce(rows):
        n, p, h = next(order)
        pivots, loss = pivot_rows(rows, sum(UNCUT[n-1]), UNCUT[n][h], jet)
        losses.append(loss)
        if p == h == 3:
            assert max(losses) == (0, 1, 1, 2, 2)[n-1]
            jet.precision -= max(losses)
            losses.clear()
        return pivots

    basis, image = field_paths(windows, jet.constant(1), reduce)
    centers = [[] for _ in range(4)]
    for h, dimension in enumerate((3, 2, 0, 0)):
        rows, width = central_equations(basis, image, h)
        free, projection, loss = saturated_projection(rows, width, dimension, jet)
        offset = 0
        for p in range(4):
            words = basis[p, 4, h]
            coordinates = tuple(col-offset for col in free if offset <= col < offset+len(words))
            if coordinates:
                normalization.append((h, p, tuple(words), coordinates))
            offset += len(words)
        assert loss == (0 if h == 0 else 1)
        for j in range(1 if h == 0 else 0, dimension):
            offset = 0
            for p in range(4):
                paths = basis[p, 4, h]
                centers[p].append((h, {word: projection[offset+i][j] for i, word in enumerate(paths)
                                       if j in projection[offset+i]}))
                offset += len(paths)
    jet.precision -= 1
    assert jet.precision == precision-7
    tables = initial_products(jet, basis, image, centers, (0, 0, 0, 0, 0, 2))
    return tables, normalization


def verify():
    s = poly([0, 1])
    assert fmpz_mod_poly_ctx(3)([1, 1, 0, 2, 1]).is_irreducible()
    assert (1+s*s)**4 + (2*s*(1-s*s))**2 == MODULUS*poly([1, 2, 0, 1, 1])
    exact = IntegralPoleJet
    exact.precision = 10
    z0, z1, z = fmpz_poly([]), fmpz_poly([1]), fmpz_poly([0, 1])
    _, exact_coordinates = initial_tables(exact, 10, exact([z1], fmpz_poly([2])),
                                          exact([z0, z1]), exact([z]))
    assert exact.precision == 1
    jet = PoleJet
    jet.precision = 17
    (dims, left, right), residue_coordinates = initial_tables(
        jet, 17, jet.constant(2), jet([ZERO, ONE]), jet([poly([0, 1])]))
    assert exact_coordinates == residue_coordinates and jet.precision == 8
    for n in range(6, 14):
        losses = [overlap(n, p, h, dims, left, right, jet) for p, h in product(range(4), repeat=2)]
        assert max(losses) == (2, 1, 1, 1, 1, 0, 1, 0)[n-6]
        jet.precision -= max(losses)
        assert jet.precision > 0
        print("Pole chart degree", n, "remaining precision", jet.precision, flush=True)


if __name__ == "__main__":
    verify()
