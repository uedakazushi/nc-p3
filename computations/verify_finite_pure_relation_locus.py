"""Check the elementary pencil statements in the finite-pure-locus proof."""
import sympy as sp


def verify():
    s, t = sp.symbols('s t')
    left = sp.Matrix([[s, 0], [t, s], [0, t]])
    kernel = sp.Matrix([[t*t, -s*t, s*s]])
    assert kernel*left == sp.zeros(1, 2)
    a = sp.symbols('a:6')
    linear = sp.Matrix([[a[i]*s+a[i+3]*t for i in range(3)]])
    equations = [c for p in linear*left for c in sp.Poly(p, s, t).coeffs()]
    assert sp.linear_eq_to_matrix(equations, a)[0].rank() == 6
    singular = sp.Matrix([[s, 0, 0], [t, 0, 0], [0, s, t]])
    assert singular.rank() == 2 and singular.det() == 0
    # The last two columns are independent pure tensors with the same row.
    assert singular[:, 1:].rank() == 1
    assert sp.Matrix([s, t]).jacobian([s, t]).rank() == 2
    print("PASS: quadratic kernel, six independent coefficients, and the fixed-factor pencil")


if __name__ == '__main__':
    verify()
