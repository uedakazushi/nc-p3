"""Check both cubic jet extensions in the (11,11,131) frame over Q."""
from sympy import Matrix, symbols, zeros
from sympy.polys.matrices import DomainMatrix


def edge(first, second, reverse=False):
    relations = []
    for family, coefficients in enumerate([first, second]):
        for j in range(3):
            relation = zeros(4)
            if family:
                relation[j, 2] = -1
                relation[j + 1, 1] = 1
            else:
                relation[j + 1, 0] = 1
            for i in range(4 - j):
                relation[j + i, 3] = coefficients[i]
            if reverse:
                relation = relation[::-1, [2, 1, 0, 3]]
            relations.append(relation.reshape(1, 16))
    return Matrix.vstack(*relations)


def overlap_dimension(forward, backward):
    # The columns index monomials in cubic, successor and predecessor.
    def coefficient(row, word):
        cubic, pair = divmod(word, 16)
        successor, predecessor = divmod(pair, 4)
        relation, fixed = divmod(row % 24, 4)
        if row < 24:
            return forward[relation, 4*cubic + successor] if fixed == predecessor else 0
        return backward[relation, 4*cubic + predecessor] if fixed == successor else 0

    system = DomainMatrix.from_Matrix(Matrix(48, 64, coefficient))
    return 48 - system.to_field().rank()


def mixed_coefficients():
    x = symbols('x')
    e0 = Matrix(2, 2, symbols('a:4'))
    e1 = Matrix(2, 2, symbols('b:4'))
    h0, h2, h3 = [Matrix(1, 2, symbols(name)) for name in ('r:2', 'u:2', 'w:2')]
    j0, j2, j3 = [Matrix(2, 1, symbols(name)) for name in ('s:2', 'v:2', 'z:2')]
    h = h0 + x**2*h2 + x**3*h3
    j = j0 + x**2*j2 + x**3*j3
    forward = (h * (e1*j2 + x*e0*j0 + x**2*e1*j0))[0].expand()
    reverse = ((h2*e0 + x*h0*e1 + x**2*h0*e0) * j)[0].expand()
    expected = [-h2*e0*j3, -h2*e0*j2, h2*e1*j2, h3*e1*j2]
    for i, product in enumerate(expected):
        assert (forward.coeff(x, i) - reverse.coeff(x, 3 - i) - product[0]).expand() == 0


def main():
    mixed_coefficients()
    cases = [
        ((2, 5), (3, 7), (1, 1), (2, 3), 0),
        ((0, 0), (1, 2), (2, 3), (1, 6), 3),
        ((1, 2), (2, 4), (5, 1), (10, 2), 4),
        ((1, 0), (0, 0), (1, 0), (0, 0), 4),
    ]
    for h2, h3, j2, j3, expected in cases:
        forward = edge([1, 0, h2[0], h3[0]], [0, 0, h2[1], h3[1]])
        backward = edge([1, 0, j2[0], j3[0]], [0, 0, j2[1], j3[1]], reverse=True)
        assert forward.rank() == backward.rank() == 6
        assert overlap_dimension(forward, backward) == expected
    print('PASS: four symbolic mixed coefficients and exact cubic overlap ranks')


if __name__ == '__main__':
    main()
