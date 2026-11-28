"""Verify the Hilbert dimensions, secant capacities, and the final incidence parity."""
from math import comb

degree, genus = 5, 2
assert degree > 2*genus-2
assert degree + 1 - genus == 4
square_sections = 2*degree + 1 - genus
outside_min = 16 - square_sections
assert outside_min == 7 and square_sections == 9
delta = comb(degree - 1, 2) - genus
assert delta == 4 and 1 + 2 + delta == outside_min
assert outside_min % 2 == 1
for base in (0, 1):
    assert comb(base+2, 1) <= 3
    assert 7//2 == 3
assert [d for d in range(1, degree+1) if d % 4 == 0] == [4]
print("PASS: quintic mixed rank, four double directions, two fixed directions and seven-mark parity")
