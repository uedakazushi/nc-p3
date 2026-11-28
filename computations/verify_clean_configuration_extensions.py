"""Check both orientations' folding obstruction on every edge of the rational configuration."""

import sympy as sp

from verify_klein_elliptic_completion import four_relations, rational_example
from verify_opposite_sextic_threshold import alternation_rank


def main():
    window, _, reverse = rational_example()
    transposition = [4*j+i for i in range(4) for j in range(4)]
    for block in four_relations(window, reverse):
        edge = sp.Matrix([coefficients for _, coefficients in block]).T
        assert alternation_rank(edge) == alternation_rank(edge[transposition, :]) == 24
    print("All eight oriented edges have invertible alternation maps and admit no symmetric folding.")


if __name__ == "__main__":
    main()
