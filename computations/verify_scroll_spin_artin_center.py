"""Compute the minimal relations of the truncated half-canonical section ring."""

WEIGHTS = (1,)*4 + (3,)*4
UNIT = (0,)*8


def monomial_value(exponents):
    nx, ny = sum(exponents[:4]), sum(exponents[4:])
    j = sum((i % 4)*n for i, n in enumerate(exponents))
    s = 3*(nx+ny) - j + nx//2
    t = j + nx//2
    if ny >= 2 or max(s, t) >= 8:
        return None
    return s, t, nx % 2, ny


def multiply(m, i):
    return m[:i] + (m[i]+1,) + m[i+1:]


def main():
    by_degree = [[UNIT]]
    relations, hilbert = [], []
    for degree in range(10):
        if degree:
            by_degree.append(sorted({multiply(m, i) for i, w in enumerate(WEIGHTS)
                                     if w <= degree for m in by_degree[degree-w]}))
        monomials = by_degree[-1]
        parent = {m:m for m in monomials}
        parent[None] = None

        def root(m):
            while parent[m] != m:
                parent[m] = parent[parent[m]]
                m = parent[m]
            return m

        # Inherit every lower-weight linear relation after multiplying by one generator.
        for i, weight in enumerate(WEIGHTS):
            if weight > degree:
                continue
            first = {None: None}
            for m in by_degree[degree-weight]:
                value, target = monomial_value(m), multiply(m, i)
                if value in first:
                    parent[root(target)] = root(first[value])
                else:
                    first[value] = target
        inherited = len({root(m) for m in monomials} - {root(None)})
        dimension = len({monomial_value(m) for m in monomials} - {None})
        relations.append(inherited - dimension)
        hilbert.append(dimension)
    assert hilbert == [1,4,7,8,7,4,1,0,0,0]
    assert relations == [0,0,3,6,10,6,10,0,0,0]
    for degrees, bounds in [((4,4,8,8), (22,21,24)), ((4,4,4,8), (18,17,20))]:
        top = sum(d-1 for d in degrees)
        assert (max(top+2, *degrees), top+1, top+4) == bounds
    print("Artin Hilbert vector, minimal relations in weights two through six, and central-module cutoffs checked.")


if __name__ == "__main__":
    main()
