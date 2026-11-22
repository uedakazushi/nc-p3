"""Rational certificates for a non-bisymmetric four-window, including its parameter test."""
from itertools import combinations, product
from math import comb

import sympy as sp


I = sp.eye(4)
J = sp.diag(-1, 1, 1, 1)
tensor = sp.kronecker_product
pairs = list(product(range(4), repeat=2))
flip = [4*j+i for i, j in pairs]


def relation_value(x, f, sign):
    u, v = x[1:, :], f[1:, :]
    mixed = sp.Matrix([u[(j+1)%3]*v[(j+2)%3] + u[(j+2)%3]*v[(j+1)%3]
                       for j in range(3)])
    return (x[0]*v + sign*u.cross(v)).col_join(f[0]*u + sign*mixed)


def gram_half(sign):
    return sp.Matrix.hstack(*(relation_value(I[:, i], I[:, j], sign) for i, j in pairs))


def verify_geometry(b, g):
    w = b.T*g
    assert all(w[4*i+j, 4*k+l] == w[4*k+j, 4*i+l] for i,j,k,l in product(range(4), repeat=4))
    opposite = sp.Matrix([[w[4*i+j, 4*k+l] for j,l in pairs] for i,k in pairs])
    assert b.rank() == g.rank() == 6 and opposite.rank() == 10
    for cut in range(4):
        one = sp.zeros(4,64)
        for i,j,k,l in product(range(4), repeat=4):
            letters = [i,j,k,l]
            row = letters.pop(cut)
            one[row,16*letters[0]+4*letters[1]+letters[2]] = w[4*i+j,4*k+l]
        assert one.rank() == 4
    r, x, y, z = sp.symbols('r x y z')
    evaluation = sum((a*b[:, 4*i:4*i+4] for i, a in enumerate((r,x,y,z))), sp.zeros(6,4))
    minors = [evaluation[list(rows), :].det() for rows in combinations(range(6), 4)]
    chart = sp.groebner([m.subs(r, 1) for m in minors], x,y,z)
    expected = [x*x-4*y*y*z*z-y*y-z*z, 2*x*y*y+2*x*z*z+x, x*y*z,
                2*y**3+y, 2*z**3+z]
    assert chart == sp.groebner(expected, x,y,z)
    assert sp.groebner([m.subs({r:0,x:1}) for m in minors], y,z).is_zero_dimensional
    print("Adjacent/opposite ranks 6/10; finite even projection with 13+7 geometric points.")


def quotient_map(relations):
    # Relations are columns. Nonpivot coordinate vectors give the quotient basis.
    reduced, pivots = relations.T.to_DM().rref()
    free = [i for i in range(relations.rows) if i not in pivots]
    result = sp.zeros(len(free), relations.rows)
    for i, f in enumerate(free):
        result[i, f] = 1
    if pivots:
        reduced = reduced.to_Matrix()
        for j, v in enumerate(pivots):
            result[:, v] = -reduced[j, free].T
    return result


def path_maps(relations, length):
    mul = {(p, 1): [I[:, i] for i in range(4)] for p in range(2)}
    dimensions = {0:[1,1], 1:[4,4]}
    for n in range(2, length+1):
        dims = []
        for p in range(2):
            tail = mul[1-p, n-1]
            rows = relations[p]
            blocks = [[sum((rows[k,4*i+j]*tail[j] for j in range(4) if rows[k,4*i+j]),
                           sp.zeros(*tail[0].shape)) for k in range(rows.rows)] for i in range(4)]
            q = quotient_map(sp.Matrix.vstack(*(sp.Matrix.hstack(*row) for row in blocks)))
            size = tail[0].rows
            mul[p,n] = [q[:,i*size:(i+1)*size] for i in range(4)]
            dims.append(q.rows)
        dimensions[n] = dims
    return mul, dimensions


def intersections(a, b):
    candidates = a.row_join(-b).nullspace()
    return a*sp.Matrix.hstack(*(x[:a.cols, :] for x in candidates))


def verify_parameters(mul):
    ev = [(0,1,1,1), (-1,0,2,-2), (-1,-2,0,2), (-1,2,-2,0)]
    od = [(0,1,1,1), (-1,0,-2,2), (-1,2,0,-2), (-1,-2,2,0)]
    assert sp.Matrix(ev).det() and sp.Matrix(od).det()
    for p in range(2):
        left, right = (ev,od) if p == 0 else (od,ev)
        for a in range(4):
            for i in range(4):
                chi = 1 if a == 0 or i in (0,a) else -1
                lhs = sum((left[a][j]*mul[p,3][j]*mul[1-p,2][j]*I[:,i]
                           for j in range(4)), sp.zeros(mul[p,3][0].rows,1))
                rhs = sum((right[a][j]*mul[p,3][i]*mul[1-p,2][j]*I[:,j]
                           for j in range(4)), sp.zeros(lhs.rows,1))
                assert lhs == chi*rhs


def verify_short_complex(relations, mul):
    spaces = {(p,1): I for p in range(2)}
    for p in range(2):
        spaces[p,2] = relations[p].T
    for n in (3,4):
        for p in range(2):
            spaces[p,n] = intersections(tensor(spaces[p,n-1],I), tensor(I,spaces[1-p,n-1]))
            assert spaces[p,n].cols == comb(4,n)
    for p in range(2):
        ranks = []
        incoming = []
        for h in range(1,5):
            # Compute in the ambient prefix tensor power, avoiding choices of coordinates in W.
            for tail_length, dest in [(2,ranks),(1,incoming)]:
                edge = (p+h-1)%2
                prepend = sp.Matrix.hstack(*mul[edge,tail_length+1])
                basis = tensor(spaces[p,h], sp.eye(prepend.cols//4))
                d = tensor(sp.eye(4**(h-1)),prepend)*basis
                dest.append(d.to_DM().rank())
        assert ranks == [20,45,36,10] and incoming == [10,20,15,4]
        print("Koszul ranks at start", p, ranks, incoming)


if __name__ == '__main__':
    b, g = gram_half(1), gram_half(-1)
    verify_geometry(b,g)
    assert g == -b*tensor(J,J)
    relations = [b, b[:,flip]*tensor(I,J)]
    mul, dims = path_maps(relations,6)
    assert all(d == [comb(n+3,3)]*2 for n,d in dims.items())
    verify_parameters(mul)
    squares = sp.eye(16)[[0,5,10,15],:]
    _, quotient_dims = path_maps([r.col_join(squares) for r in relations],5)
    assert all(d == [comb(4,n)]*2 for n,d in quotient_dims.items())
    print("Path dimensions:", dims, "square-zero quotient:", quotient_dims)
    verify_short_complex(relations,mul)
