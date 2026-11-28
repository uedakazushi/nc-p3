"""Hilbert-scheme integrals on the two plane blowup models of a cubic kernel."""
from verify_cofactor_nineteen_clusters import intersection


def plane_blowups(degree, multiplicity, centers, weights):
    points = [(*(w - v for w in weights if w != v), -degree * v) for v in weights]
    for _ in range(centers):
        u, v, h = points.pop(0)
        points.extend([(u, v - u, h - multiplicity * u),
                       (u - v, v, h - multiplicity * v)])
    return points


for weights in [(0, 17, 271), (0, 19, 307)]:
    for degree, multiplicity, centers, expected in [(4, 1, 8, 13), (6, 2, 7, 14)]:
        assert degree**2 - centers * multiplicity**2 == 8
        assert centers * multiplicity - 3 * degree == -4
        assert intersection(plane_blowups(degree, multiplicity, centers, weights)) == expected
print("PASS: plane blowups give isolated-triple bounds 13 and 14")
