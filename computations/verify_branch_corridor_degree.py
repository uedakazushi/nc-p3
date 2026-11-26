"""Check the telescoping and boundary terms of the curve-degree formula symbolically."""
from sympy import expand, symbols

for m in range(2, 15):
    d = dict.fromkeys(range(-1, m+3), 0)
    d[1] = d[m] = 1
    d.update(zip(range(2, m), symbols(f"d2:{m}")))
    differences = [d[i+2]-3*d[i+1]+3*d[i]-d[i-1] for i in range(m+1)]
    assert expand(sum(differences)) == 0
    assert expand(sum(i*diff for i, diff in enumerate(differences))) == m+2
print("PASS: unweighted and index-weighted degree identities for symbolic interiors")
