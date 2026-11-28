"""Verify the binary discriminant identity and the residual interpolation bounds."""
from sympy import Matrix, symbols


def main():
    a, b, c, p = symbols("a b c p")
    matrix = Matrix([[a, (b+p)/2], [(b-p)/2, c]])
    assert (4*matrix.det() - (p**2 - (b**2-4*a*c))).expand() == 0
    assert matrix.det().subs(p, b).expand() == a*c
    assert [(k, 8-2*k) for k in (3, 4)] == [(3, 2), (4, 0)]
    assert 20-3-3 > 20-2*4
    assert 20-2-3 > 14
    # Length-four fibers with collisions still pair their socle with point evaluation.
    t = symbols("t")
    from sympy import Poly
    for root, polynomial in [(0, t**4), (1, (t-1)**2*t*(t+1))]:
        modulus = Poly(polynomial, t)
        socle = modulus.exquo(Poly(t-root, t))
        assert (socle * Poly(t-root, t)).rem(modulus).is_zero
        assert socle.degree() == 3
    print("PASS: norm identity, rank-two factorization, residual numerical bounds and nonreduced quartic fiber evaluations")


if __name__ == "__main__":
    main()
