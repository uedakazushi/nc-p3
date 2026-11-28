"""Incidence codimensions, graph degree and the residual point budget."""
from math import comb

# W has dimension four in six-space. Compare dim(W intersect ker A) to j.
def incidence_codimension(kernel_dimension, j):
    return j*(6-4-kernel_dimension+j)

assert [incidence_codimension(k,1) for k in (1,2,3)] == [2,1,0]
assert incidence_codimension(3,2) == 2
assert incidence_codimension(2,2) == 4
assert comb(6,3) == 20
for curve_degree in (1, 2):
    base_curve_degree = 2*curve_degree+1
    for outside_plane in (18, 19, 20):
        other_marks = outside_plane-3*curve_degree
        residual = 20-base_curve_degree-1-other_marks
        assert other_marks > 0 and residual <= 20-outside_plane <= 2
print("PASS: Grassmann incidence and residual isolated-point bounds for both curve degrees")
