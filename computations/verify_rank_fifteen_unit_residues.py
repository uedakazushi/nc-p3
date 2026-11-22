"""Integral presentations at all four pairs (x,y) in the three-adic unit residue grid."""

from flint import fmpz_mod_poly_ctx, fq_default_ctx

from verify_rank_fifteen_cubic_zero_disk import tensor_relations
from verify_rank_fifteen_parameter_orbits import pattern
from verify_rank_fifteen_pencil_linear_residues import graded_profile
from verify_rank_fourteen_geometric_reduction import Scalar

# x, y, modulus, and K,H,n,P_+,P_-,Q_+,Q_- in the polynomial basis (base three).
POINTS = (
    (1, 1, (2, 2, 1, 0, 2, 0, 1), (459, 189, 616, 469, 649, 237, 596)),
    (1, 2, (2, 0, 0, 2, 1), (73, 1, 4, 70, 11, 11, 36)),
    (2, 1, (2, 0, 0, 2, 1), (73, 74, 33, 73, 0, 31, 74)),
    (2, 2, (2, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1),
     (116398, 507099, 150968, 503407, 439618, 161817, 85025)),
)


def element(field, encoded):
    coefficients = []
    while encoded:
        encoded, coefficient = divmod(encoded, 3)
        coefficients.append(coefficient)
    return field(coefficients)


def verify():
    for x, y, modulus, coordinates in POINTS:
        m = fmpz_mod_poly_ctx(3)(list(modulus))
        assert m.is_irreducible()
        field = fq_default_ctx(modulus=m, var="a")
        k, h, n, p, u, q, v = (element(field, c) for c in coordinates)
        U, V = k*k, h*h
        assert U != 1 and V*(U+3)/(U-1)**2 == field(x) and V**3/(U-1)**4 == field(y)
        assert n**4+(1-U)*n*n+V == field.zero() and not n.is_zero()
        assert (p*p, u*u, q*q, v*v) == (2*(k+1+h), 2*(k+1-h), 2*(k-1+h), 2*(k-1-h))
        gamma = pattern(field.one(), (k-1)/2, n/2, h/(2*n),
                        (u-p)/4, (u+p)/4, (v-q)/4, (v+q)/4)
        entries = tuple(Scalar(field.zero() if c == 0 else c) for c in gamma)
        windows = tensor_relations(entries, Scalar(field.zero()))
        choices = ((1, 0, 0), (0, 1, 0)) if (x, y) == (1, 2) else ((0, 1, 0), (0, 0, 1))
        assert graded_profile(windows, Scalar(field.one()), choices)
        print("Verified unit residues", (x, y), "in extension degree", m.degree(), flush=True)


if __name__ == "__main__":
    verify()
