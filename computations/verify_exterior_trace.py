#!/usr/bin/env python3
"""Exact four-window trace-matrix test for the exterior Koszul dual of k[x1,...,x4].

Checks the sign convention P_i^T=Q_i for the Frobenius-normalized
four-periodic marking. This is a worked example, not a verification
of AS-regularity for arbitrary quadratic relations.
"""
from fractions import Fraction
import json

DEG = 4
FULL = (1 << DEG) - 1

def wedge(a: int, b: int):
    if a & b:
        return 0, 0
    crossings = sum(1 for u in range(DEG) for v in range(DEG)
                    if (a & (1 << u)) and (b & (1 << v)) and u > v)
    return (-1 if crossings % 2 else 1), a | b

def matrix(i, trace):
    # u_a=e_a; v_b=(-1)^b e_0...hat(e_b)...e_3, so e_b*v_b=omega.
    p = [[Fraction(0) for _ in range(DEG)] for _ in range(DEG)]
    q = [[Fraction(0) for _ in range(DEG)] for _ in range(DEG)]
    for a in range(DEG):
        for b in range(DEG):
            v = FULL ^ (1 << b)
            orientation = (-1)**b
            sign, monomial = wedge(1 << a, v)
            if monomial == FULL:
                p[a][b] = trace[i % 4] * orientation * sign
            sign, monomial = wedge(v, 1 << a)
            if monomial == FULL:
                q[b][a] = trace[(i + 1) % 4] * orientation * sign
    return p, q

def transpose(a):
    return [list(x) for x in zip(*a)]

def verify():
    const = [Fraction(1)] * 4
    alternating = [Fraction((-1)**i) for i in range(4)]
    unnormalized = []
    normalized = []
    for i in range(4):
        p0, q0 = matrix(i, const)
        p1, q1 = matrix(i, alternating)
        assert p0 == [[int(r == c) for c in range(DEG)] for r in range(DEG)]
        assert q0 == [[-int(r == c) for c in range(DEG)] for r in range(DEG)]
        assert transpose(p1) == q1, f"canonical trace identity fails at {i}"
        unnormalized.append(-1)
        normalized.append(1)
    t = [Fraction(2), Fraction(3), Fraction(5), Fraction(7)]
    altered_eigenvalues = []
    for i in range(4):
        p, q = matrix(i, t)
        assert p[i % DEG]
        eigen = -t[i] / t[(i + 1) % 4]
        assert all(p[a][b] == (eigen * q[b][a]) for a in range(DEG) for b in range(DEG))
        altered_eigenvalues.append(str(eigen))
    product = Fraction(1)
    for e in altered_eigenvalues:
        product *= Fraction(e)
    assert product == 1
    return {
        "result": "PASS",
        "example": "exterior Koszul dual of commutative polynomial 4-space",
        "vertex_cases": 4,
        "pairing_matrix_size": [4, 4],
        "constant_trace_nakayama": unnormalized,
        "alternating_trace_nakayama": normalized,
        "arbitrary_rescaling_eigenvalues": altered_eigenvalues,
        "four_step_scalar_product": str(product),
        "scope": "exact example and trace-scaling identity only",
    }

if __name__ == "__main__":
    print(json.dumps(verify(), ensure_ascii=False, indent=2))
