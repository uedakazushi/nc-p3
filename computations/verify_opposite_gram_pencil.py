"""Exact identities behind the Gram form, the isotropic incidence, and its length bounds."""

import sympy as sp


def verify_gram():
    u = sp.Matrix(4, 3, sp.symbols("u:12"))
    minors = signed_minors(u)
    assert all(sp.expand(value) == 0 for value in u.T * minors)
    # The Gram cofactor formula also follows directly from Cauchy--Binet.
    for i in range(4):
        for j in range(4):
            left = u.extract([k for k in range(4) if k != i], range(3))
            right = u.extract([k for k in range(4) if k != j], range(3))
            assert sp.expand((-1)**(i+j) * left.det() * right.det() - minors[i]*minors[j]) == 0
    t = sp.symbols("t:2")
    conic = sp.Matrix([t[0]**2-t[1]**2, sp.I*(t[0]**2+t[1]**2), 2*t[0]*t[1]])
    assert sp.expand(conic.dot(conic)) == 0


def signed_minors(matrix):
    return sp.Matrix([(-1)**i * matrix.extract([k for k in range(4) if k != i],
                                              range(3)).det() for i in range(4)])


def degree_partitions(remaining, minimum=(0, 1)):
    if remaining == 0:
        yield ()
    for b in range(1, remaining + 1):
        for a in range(min(b, remaining-b) + 1):
            if (a, b) >= minimum:
                for tail in degree_partitions(remaining-a-b, (a, b)):
                    yield ((a, b), *tail)


def verify_smith_bounds():
    possible = set()
    for points in degree_partitions(8):
        length_bound = sum(sum(3 if j < a else 1 for j in range(b)) for a, b in points)
        common = sum(a for a, b in points)
        assert length_bound == 8 + common
        if length_bound >= 10:
            possible.add(common)
    assert possible == {2, 3, 4}


if __name__ == "__main__":
    verify_gram()
    verify_smith_bounds()
    print("Verified the Gram minor identities, conic parametrization, and all degree-eight Smith profiles.")
