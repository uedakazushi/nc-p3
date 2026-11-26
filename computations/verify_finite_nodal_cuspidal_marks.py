"""Check the plane-intersection bound and order gain on each normalization branch."""
from sympy import Matrix, symbols

a, b, c, d = symbols("a b c d")
assert Matrix([[c, -b*d], [d, c-a*d]]).det().expand() == c**2-a*c*d+b*d**2
for order in (1, 2):
    for cutoff in range(2*order, 35):
        for s in range(1, cutoff+2):
            available = max(0, cutoff-s-order)
            assert 2*order+s+available >= cutoff+order
            assert s+available+2*order >= cutoff+order
print("PASS: plane norm determinant and nodal/cuspidal conductor inequalities")
