"""Numerics of the Clifford equality case and intersection bounds on the kernel graph."""
from math import comb
assert comb(7,3) == 35 == 1+2*17
genus = 35
excess = 80-6*17-(1-genus)-(6-3)
assert excess == 9
assert (20-excess)*(2*4-4) > genus
assert [(d,24-d,(24-d)**2//8) for d in range(18,22)] == [
       (18,6,4),(19,5,3),(20,4,2),(21,3,1)]
print("PASS: genus and local normalization defect, sectional degrees and Hodge caps")
