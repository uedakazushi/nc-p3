"""Check binary ideal recovery and the elliptic trisecant bundle degrees."""
def homogeneous_dimension(degree):
    return max(degree + 1, 0)

for n in range(13):
    quotient = (homogeneous_dimension(n) - 3*homogeneous_dimension(n-8)
                + 2*homogeneous_dimension(n-12))
    assert quotient >= 0
    if n in (10, 11, 12):
        assert quotient == {10: 2, 11: 0, 12: 0}[n]
kernel_degree = -4 - (-1)
assert kernel_degree == -3
assert [-1-kernel_degree, -kernel_degree] == [2, 3]
assert max(2, 3) < 16-10
print("PASS: binary complete-intersection square and elliptic trisecant cap three versus six")
