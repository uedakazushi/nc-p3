"""Symbolic checks for the residual quadric and the order-three rank obstruction."""
from itertools import product

from sympy import Matrix, expand, factor, symbols, zeros


def quadric(b):
    return (b[0, 0]*b[2, 2]-(b[0, 1]+b[1, 0])*(b[1, 2]+b[2, 1])+
            b[1, 1]*(b[0, 2]+b[2, 0]+b[1, 1]))


def main():
    x, y = symbols("x:4"), symbols("y:4")
    b = Matrix(3, 3, lambda i, j: x[i]*y[j+1]-x[i+1]*y[j])
    assert expand(b.det()) == expand(quadric(b)) == 0
    p, r, q, s, h, k = symbols("p r q s h k")
    def point(t):
        return Matrix([1, t, t*t])
    test = h*point(p)*point(q).T+k*point(r)*point(s).T
    assert expand(quadric(test)-h*k*(p-r)*(p-s)*(q-r)*(q-s)) == 0

    a, b, c, d, u = symbols("A B C D u")
    terms = {(0, 0, 0, 0): 1}
    terms.update({tuple(int(j != i) for j in range(4)): z for i, z in enumerate([a, b, c, d])})
    pairs = list(product(range(3), repeat=2))
    cut = zeros(9)
    for bits, coeff in terms.items():
        for i, sign_i in [(2, 1), (1, -1)]:
            for j, sign_j in [(0, u), (3, -1)]:
                index = list(bits)
                index[i] += 1
                index[j] += 1
                cut[pairs.index(tuple(index[:2])), pairs.index(tuple(index[2:]))] += coeff*sign_i*sign_j
    blocks = [cut.extract([i for i, z in enumerate(pairs) if sum(z) % 3 == w],
                          [i for i, z in enumerate(pairs) if (sum(z)+w) % 3 == 2]) for w in range(3)]
    assert [factor(m.det()) for m in blocks] == [c*d*u, -u*(a*c*u+b*d), a*b*u**2]
    assert blocks[0] == Matrix([[0, -1, 0], [c, d-a*u, 0], [0, c*u-b*u, d*u]])
    assert blocks[2] == Matrix([[0, a, 0], [-u, b-c, a*u-d], [0, 0, b*u]])
    print("PASS: residual quadric and determinant, secant fiber, and all three weight blocks")


if __name__ == "__main__":
    main()
