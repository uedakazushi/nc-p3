"""Reducedness and normality certificate on the irreducible common-quadric locus."""

from math import prod

import sympy as sp
from flint import nmod_mat, nmod_mpoly_ctx

from verify_normal_configuration_octics import smooth_plane_form
from verify_opposite_root_triple_octics import monomials, octic_evaluation, quotient_projection
from verify_semisimple_period_obstruction import P, inverse
from verify_single_edge_instanton_extension import coefficients_to_field
from verify_trisecant_middle_simple_obstruction import common_quadric_transform, middle_relations


def main():
    relation = sp.Matrix(coefficients_to_field(middle_relations(common_quadric_transform())).tolist())
    quartics, projection = quotient_projection(relation, 4)

    def multiply(i):
        return projection[:, [quartics[tuple(sorted((i,)+word))] for word in monomials(3)]]

    affine = (inverse(multiply(0))@multiply(1)) % P
    characteristic = nmod_mat(affine.tolist(), P).charpoly()
    assert characteristic.degree() == 20 and characteristic.gcd(characteristic.derivative()).degree() == 0
    equations = nmod_mat(octic_evaluation(relation).tolist(), P)
    basis, nullity = equations.nullspace()
    assert nullity == 1
    ctx = nmod_mpoly_ctx.get(("x", "y", "z"), P)
    x, y, z = ctx.gens()
    coordinates = (x, y, z, x+2*y+5*z)
    section = sum(int(basis[i, 0])*prod(coordinates[j] for j in word)
                  for i, word in enumerate(monomials(8)))
    smooth_plane_form(section, P)
    print("The twenty middle points are reduced, with a unique triple octic and a smooth plane section.")


if __name__ == "__main__":
    main()
