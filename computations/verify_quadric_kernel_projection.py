"""Recover the symmetric root of a rational (2,2) projection by exact linear algebra."""
from itertools import combinations_with_replacement
from sympy import Matrix

QUADRICS = [
    [-11,10,11,-11,-13,-3,1,-18,-18],
    [16,-14,12,17,9,11,-16,15,-18],
    [-19,-4,10,-15,1,7,-7,13,5],
    [9,2,7,-8,15,3,-4,3,-8],
]
PAIRS = list(combinations_with_replacement(range(4), 2))


def monomials(a, b):
    return [(i,j) for i in range(a+1) for j in range(b+1)]


def basis(a, b):
    return [{key: 1} for key in monomials(a,b)]


def multiply(left, right):
    result = {}
    for (i,j),v in left.items():
        for (k,l),w in right.items():
            key = (i+k,j+l)
            result[key] = result.get(key,0)+v*w
    return result


def evaluations(polys, a, b):
    return Matrix([[f.get(key, 0) for f in polys] for key in monomials(a,b)])


def forms():
    return [dict(zip(monomials(2,2), entries)) for entries in QUADRICS]


def symmetric_root():
    v, y = forms(), basis(1,1)
    linear = evaluations([multiply(h,l) for h in v for l in y],3,3)
    assert linear.to_DM().rank() == 16
    products = [multiply(multiply(v[i],v[j]),l) for i,j in PAIRS for l in y]
    relations = evaluations(products,5,5).to_DM().to_field().nullspace().to_Matrix()
    assert relations.shape == (4,40)
    coeffs = [relations[:,4*k:4*k+4].T for k in range(10)]
    constraints = Matrix([
        [c[i,k]*int(l==j)-c[j,k]*int(l==i) for k in range(4) for l in range(4)]
        for c in coeffs for i in range(4) for j in range(i)
    ])
    solution = constraints.to_DM().to_field().nullspace().to_Matrix()
    assert solution.shape == (1,16)
    change = solution.reshape(4,4)
    assert change.det() != 0
    coeffs = [c*change for c in coeffs]
    assert all(c.is_symmetric() for c in coeffs)
    entries = Matrix([[c[i,j] for i,j in PAIRS] for c in coeffs])
    assert entries.to_DM().rank() == 10
    x = [sum(h.values()) for h in v]
    value = sum((x[i]*x[j]*c for (i,j),c in zip(PAIRS,coeffs)), Matrix.zeros(4))
    assert value.rank() == 3 and value*Matrix.ones(4,1) == Matrix.zeros(4,1)
    return entries


if __name__ == "__main__":
    symmetric_root()
    print("PASS over Q: four quadratic syzygies, unique symmetric presentation, full span, birational projection")
