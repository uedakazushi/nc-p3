"""All degree-two and degree-three unit residues on the vertical x=1 chart."""

from itertools import product

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx

from verify_rank_fifteen_unit_residues import element, verify_unit_point

# Polynomials are in ascending degree. Coordinates use the base-three encoding
# of the coefficients in the indicated extension's polynomial basis.
FIELDS = {
    6: (2, 2, 1, 0, 2, 0, 1),
    8: (2, 2, 2, 0, 1, 2, 0, 0, 1),
    12: (2, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1),
    18: (2, 0, 2, 0, 2, 1, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1),
}

# Minimal polynomial of y, field degree, y, then K,H,n,P+,P-,Q+,Q-.
FIBERS = (
    ((1, 0, 1), 12, 483498, (435704, 344254, 338756, 83179, 494561, 57734, 253525)),
    ((2, 1, 1), 8, 3166, (4797, 5779, 3601, 5732, 48, 6170, 1953)),
    ((2, 2, 1), 8, 6323, (4798, 616, 5105, 1009, 3042, 1401, 5711)),
    ((1, 0, 2, 1), 12, 8029, (234298, 75446, 307314, 370208, 235788, 449572, 169545)),
    ((1, 1, 2, 1), 12, 110583, (304113, 512139, 156524, 485371, 298780, 487629, 79691)),
    ((1, 2, 0, 1), 6, 649, (614, 646, 58, 718, 425, 334, 696)),
    ((1, 2, 1, 1), 6, 135, (149, 182, 189, 425, 491, 696, 294)),
    ((2, 0, 1, 1), 18, 370171635, (305037659, 5246843, 243929955, 240904061, 375685257, 225306306, 271917414)),
    ((2, 1, 1, 1), 6, 460, (540, 596, 317, 317, 542, 469, 332)),
    ((2, 2, 0, 1), 18, 97828222, (280018502, 16383589, 173100917, 50557231, 22466366, 43226751, 318862864)),
    ((2, 2, 2, 1), 18, 256678021, (252372802, 78082063, 215643424, 195340852, 132610298, 3705705, 91907813)),
)


def verify():
    ring = fmpz_mod_poly_ctx(3)
    expected = {(*a, 1) for degree in (2, 3) for a in product(range(3), repeat=degree)
                if ring([*a, 1]).is_irreducible()}
    assert len(expected) == len(FIBERS) == 11 and expected == {m for m, _, _, _ in FIBERS}
    for m, degree, encoded_y, coordinates in FIBERS:
        modulus = ring(list(FIELDS[degree]))
        assert modulus.is_irreducible()
        field = fq_default_ctx(modulus=modulus, var="a")
        y = element(field, encoded_y)
        assert fq_default_poly_ctx(field)(list(m))(y).is_zero()
        k, h = (element(field, code) for code in coordinates[:2])
        v = 1/k
        assert (v-v**3)**2 == y and k-v == h
        verify_unit_point(field, field.one(), y, coordinates, ((0, 1, 0), (0, 0, 1)))
        print("Verified unit residue polynomial", m, "in degree", degree, flush=True)


if __name__ == "__main__":
    verify()
