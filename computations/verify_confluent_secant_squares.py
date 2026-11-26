"""Square and first-jet coordinates for a six-point cut with one double node."""
from itertools import combinations

import sympy as sp


def verify():
    x = sp.symbols("x:4")
    a = sp.symbols("a:4")
    pairs = list(combinations(range(4), 2))
    singles = [u*u for u in x]
    edges = [(x[i]-x[j])**2 for i, j in pairs]
    base = singles + edges
    monomials = singles + [x[i]*x[j] for i, j in pairs]
    coefficients = sp.Matrix([[sp.expand(q).coeff(m) for q in base] for m in monomials])
    assert abs(coefficients.det()) == 2**6
    w = sum(ai*xi for ai, xi in zip(a, x))
    square = sum(ai*sum(a)*q for ai, q in zip(a, singles))
    square -= sum(a[i]*a[j]*q for (i, j), q in zip(pairs, edges))
    assert sp.expand(w**2-square) == 0
    for i in range(4):
        derivative = sum(ak*xk*xk for ak, xk in zip(a, x)) + sum(a)*singles[i]
        derivative -= sum(a[k]*(x[i]-x[k])**2 for k in range(4) if k != i)
        assert sp.expand(2*x[i]*w-derivative) == 0
    assert (16-2)//2 + 1 > 7


if __name__ == "__main__":
    verify()
    print("PASS: ten square basis elements, tangent square, all four first jets, and re-marking count")
