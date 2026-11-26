"""Auxiliary symbolic identities in the genus-two pencil and the cubic-cone cases."""
from sympy import Matrix, expand, symbols


def main():
    a, b, c, d, s, t = symbols("a b c d s t")
    S, T = a*s+b*t, c*s+d*t
    v = Matrix([s*s, s*t, t*t])
    w = Matrix([S*S, S*T, T*T])
    delta = s*T-t*S
    # The three cofactors have precisely the fixed divisor of a nonidentity
    # PGL2 map. The remaining three coordinates are base-point-free for det != 0,
    # as v and w are the Veronese squares of two nonzero vectors.
    assert (v.cross(w)-delta*Matrix([t*T, -s*T-t*S, s*S])).applyfunc(expand) == Matrix.zeros(3, 1)
    # On F_2 the class m*section+5*fiber has the indicated arithmetic genus.
    for m, genus, contact in [(1, 0, 3), (2, 2, 1)]:
        self_intersection = -2*m*m+10*m
        canonical_intersection = -10  # K=-2*section-4*fiber
        assert 1+(self_intersection+canonical_intersection)//2 == genus
        assert 5-2*m == contact
    g = 2
    trace_degree = g+1
    assert 2*trace_degree == 6 and 2*g-2 == 2
    print("PASS: conic cofactors including double fixed roots, cone adjunction and double-cover lengths")


if __name__ == "__main__":
    main()
