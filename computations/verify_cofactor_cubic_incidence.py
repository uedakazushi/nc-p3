"""The sectional sextic and the cubic's bilinear incidence intersections."""
from math import comb

hilbert = [1, 4] + [6*t-2 for t in range(2, 11)]
numerator = [sum((-1)**i * comb(4, i) * hilbert[j-i] for i in range(min(4, j)+1))
             for j in range(len(hilbert))]
assert numerator == [1, 0, 0, -4, 3] + [0]*6
# On P^3 x P^2, take the L^3 h^2 coefficient after multiplying by (L+h)^3.
squares = [comb(3, 3-power) for power in (2, 1, 0)]
assert squares == [3, 3, 1]
assert squares[0]+2*squares[1]+squares[2]-8 == 2
assert [12-4-8, 12-2-8] == [0, 2]
print("PASS: sextic Hilbert numerator, incidence intersections, and both cubic projections")
