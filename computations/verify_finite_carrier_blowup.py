"""Check the integral blowup lattices and conductor inclusions over Q[[x]]."""
from sympy import Matrix, Symbol, diag, eye, zeros


def verify(d, coefficients):
    x = Symbol("x")
    # Columns express multiplication by z on the strict transform.
    z = zeros(d)
    for j in range(d - 1):
        z[j + 1, j] = 1
    z[:, d - 1] = Matrix([-v for v in reversed(coefficients)])
    lattice = diag(*(x**j for j in range(d)))
    y = (lattice.inv() * (x*z) * lattice).applyfunc(lambda a: a.cancel())
    assert all(a.is_polynomial(x) for a in y)
    poly = y**d
    for i, a in enumerate(coefficients, 1):
        poly += x**i * a * y**(d-i)
    assert poly.applyfunc(lambda a: a.expand()) == zeros(d)
    conducted = lattice.inv() * x**(d-1) * eye(d)
    assert all(a.is_polynomial(x) for a in conducted)
    assert not (conducted[d-1, d-1]/x).is_polynomial(x)
    # On the original lattice g=x^6 kills the local ring. Every generator
    # of m^2 kills the whole discrepancy g C' modulo g C for d <= 3.
    for i in range(3):
        product = lattice.inv() * x**2 * z**i
        assert all(a.expand().is_polynomial(x) for a in product)
    if d == 3:
        # The loss of one x before using the lower row would matter.
        assert not (lattice.inv() * x)[2, 2].is_polynomial(x)


if __name__ == "__main__":
    x = Symbol("x")
    for d in (2, 3):
        verify(d, [Symbol(f"b{i}") for i in range(d)])
        verify(d, [0]*d)
        verify(d, [0]*(d-1) + [x**7])
    print("PASS: double/triple plane carriers, lattice inclusions and the one-power discrepancy")
