"""Check the quartic line incidence graph and multiplication in the two polygon models."""
from itertools import combinations, product

from sympy import Matrix, binomial, symbols


def curves_of_lines():
    vertices = range(4)
    edges = list(combinations(vertices, 2))
    remaining = []
    for flags in product([0, 1], repeat=len(edges)):
        graph = {edge for edge, chosen in zip(edges, flags) if chosen}
        if any(sum(i in edge for edge in graph) < 2 for i in vertices):
            continue
        if any(set(combinations(face, 2)) <= graph for face in combinations(vertices, 3)):
            continue
        remaining.append(graph)
    assert len(remaining) == 3
    assert all(len(graph) == 4 for graph in remaining)


def sections(m, d, glue):
    # Coefficient of z^j on component a is recorded at (a,j).
    basis = []
    for a in range(m):
        basis.append({(a, 0): 1, ((a-1) % m, d): glue if a == 0 else 1})
        basis.extend({(a, j): 1} for j in range(1, d))
    return basis


def polygon_products(m, gluing):
    d = 4//m
    left, right = (sections(m, d, g) for g in gluing)
    columns = []
    for v, w in product(left, right):
        column = [0]*(m*(2*d+1))
        for (a, i), vi in v.items():
            for (b, j), wj in w.items():
                if a == b:
                    column[a*(2*d+1)+i+j] += vi*wj
        columns.append(column)
    return Matrix(columns).T


def main():
    curves_of_lines()
    u, v = symbols("u v", nonzero=True)
    for m in (2, 4):
        for left, right in [(u, v), (1, 1), (1, -1)]:
            assert polygon_products(m, (left, right)).rank() == 8
        assert binomial(6, 3)-4*4 == 4
    print("PASS: four-line incidence enumeration and cyclic quartic multiplication matrices")


if __name__ == "__main__":
    main()
