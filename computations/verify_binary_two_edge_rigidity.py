"""Check the compressed cubic intersection for distinct binary marks and their Artin jets."""
from functools import reduce
from itertools import combinations

import sympy as sp

from verify_binary_secant_obstructions import rank


def vector(x):
    return sp.Matrix([1, x, x**2, x**3])


def mapped_curve(xs, ys):
    left = sp.Matrix.hstack(*(vector(i) for i in xs[:4]))
    right = sp.Matrix.hstack(*(vector(i) for i in ys[:4]))
    if len(xs) == 5:
        a, b = left.inv()*vector(xs[4]), right.inv()*vector(ys[4])
        scale = [v/u for u,v in zip(a,b)]
    else:
        scale = [1, 3, 2, 7]
    return right*sp.diag(*scale)*left.inv()


def check(matrix, m):
    x = sp.symbols("x")
    a = matrix*vector(x)
    # The degree accounts also for a possible common zero at infinity.
    minors = [sp.Poly(a[i]*a[j+1]-a[i+1]*a[j], x) for i,j in combinations(range(3),2)]
    common = reduce(sp.gcd, minors)
    assert common.degree()+6-max(p.degree() for p in minors) == m
    cols = []
    for i in range(3):
        for j in range(3):
            for k in range(2):
                v = [sp.S.Zero]*3
                v[k] = -x**j*a[i+1]
                v[k+1] = x**j*a[i]
                cols.append(sp.Matrix([sp.expand(q).coeff(x,l) for q in v for l in range(6)]))
    assert 18-rank(sp.Matrix.hstack(*cols)) == m


def main():
    check(mapped_curve([0,1,2,3], [1,3,2,7]), 4)
    check(mapped_curve([0,1,2,3,4], [1,3,2,7,5]), 5)
    x = sp.symbols("x")
    v = vector(x+x**2)
    c = sp.symbols("c:5")
    product = (sum(a*x**i for i,a in enumerate(c)) * v).applyfunc(sp.expand)
    truncation = sp.Matrix([f.coeff(x,4) for f in product]).jacobian(c).nullspace()[0]
    truncation /= truncation[0]
    jet = product.subs(dict(zip(c,truncation)))
    matrix = sp.Matrix(4,4,lambda i,j:jet[i].coeff(x,j))
    assert matrix.det() != 0
    check(matrix, 5)
    print("Compressed cubic dimensions 4, 5 and 5 verified, including a single length-five intersection.")


if __name__ == '__main__':
    main()
