"""Sharp numerical inconsistency of quartic self-linkage and the octic factor bound."""
from math import comb
dimension_four = comb(7,3)
bound_genus = 4*16+1-(dimension_four-1)
assert bound_genus == 31
assert 6*16+2-2*bound_genus+5 == 41 > 2*20
for degree in range(1,22):
    least_genus = (6*degree+2+5-40+1)//2
    assert least_genus == 3*degree-16
    if dimension_four <= 4*degree+1-least_genus:
        assert degree >= 18
    if dimension_four <= 2*degree+1:
        assert degree >= 17
assert [(e,d) for e in range(5,9) for d in range(17,23)
        if d <= comb(e-1,2)] == [(8,d) for d in range(17,22)]
print("PASS: quartic self-linkage contradiction, Clifford bounds and integral octic degrees")
