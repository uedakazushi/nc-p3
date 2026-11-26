"""Check the secant determinant, all-length Euler identities, and a basepoint pencil."""

from sympy import Matrix, Poly, factor, simplify, symbols

a, b, c, d, z, t = symbols("a b c d z t")
edge = Poly((t - z) * (a + b*z + c*t + d*z*t), z, t)
matrix = Matrix(3, 3, lambda i, j: edge.nth(i, j))
assert factor(matrix.det() - (c - b)*(a*d - b*c)) == 0

# Removing endpoint brackets amounts to injecting binary quadratics into cubics.
p, q = symbols("p q")
def linear_factor(root):
    return Matrix(4, 3, lambda i, j: int(i == j + 1) - root*int(i == j))
left, right = linear_factor(p), linear_factor(q)
assert (left * matrix * right.T).extract([1, 2, 3], [1, 2, 3]).det().expand() == matrix.det().expand()

n = symbols("n", integer=True, positive=True)
def four_letter(length):
    return (length + 1)*(length + 2)*(length + 3)/6
assert simplify(four_letter(n) - 2*four_letter(n-1) + four_letter(n-2) - (n+1)) == 0
assert simplify((1/(1-t)**2 - t/(1-t)**2) - 1/(1-t)) == 0

# Ordinary polynomial window. Take the rank-two pencil x1,x3 on the
# twisted cubic (1,z,z^2,z^3), whose moving part has degree two.
# At p=0 the extra generator of the point ideal is x2. In the quotient
# multiplication by x2 sends the next line module to the previous one.
for length in range(10):
    line_basis = {(length - j, j) for j in range(length + 1)}  # exponents of x0,x2
    point_basis = {(length, 0)}
    previous_image = {(length-1-j, j+1) for j in range(length)}
    assert line_basis == point_basis | previous_image and point_basis.isdisjoint(previous_image)
    assert len(line_basis) == length + 1

# Two independent coefficients of each relation and the projected cubic syzygy.
x0, x1, x2, x3 = xs = symbols("x0:4")
def wedge(i, j):
    return Matrix(4, 1, lambda k, _: int(k == i)*xs[j] - int(k == j)*xs[i])
r, s1, s2 = wedge(1, 3), wedge(1, 2), wedge(2, 3)
assert r[2] == 0 and (s1[2], s2[2]) == (-x1, x3)
assert r*x2 - s1*x3 - s2*x1 == Matrix.zeros(4, 1)
print("PASS: secant determinant and endpoint injections, Hilbert identities, point filtration and lifted syzygy")
