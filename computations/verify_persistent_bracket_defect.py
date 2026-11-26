"""Exact nine-word determinant on the diagonal chart and the persistent cubic."""
import sympy as sp


def comm(a, b):
    return a*b-b*a


def words(a, b):
    s = a*b+b*a
    return [sp.eye(3), a, b, a*a, b*b, s, a*b*a, b*a*b, s*s]


def verify():
    x, y = sp.symbols("x y")
    a = sp.diag(x, y, -x-y)
    b = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    matrix = sp.Matrix.hstack(*(w.reshape(9, 1) for w in words(a, b)))
    assert sp.expand(matrix.det()+comm(a, b).det()**3) == 0
    # A nonsemisimple first generator.
    a = sp.Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    b = sp.Matrix([[0, 1, 2], [3, 2, -1], [1, 0, -2]])
    assert comm(a, b).det() != 0
    assert sp.Matrix.hstack(*(w.reshape(9, 1) for w in words(a, b))).det() == -comm(a, b).det()**3
    s = sp.Matrix(2, 2, sp.symbols("s:4"))
    def odd(i):
        w = sp.eye(2)[:, i]
        return sp.zeros(2).row_join(w).col_join((w.T*s).row_join(sp.zeros(1)))
    p = sp.diag(0, 0, 1)
    u, v = odd(0), odd(1)
    h = sp.Symbol("h")
    z = comm(u, v)
    assert sp.expand((h*z+x*comm(p, u)+y*comm(p, v)).det()-h**3*z.det()) == 0
    column = sp.Matrix([x, y])
    compressed = x*comm(p, u)+y*comm(p, v)
    assert compressed[:2, 2] == -column and compressed[2, :2] == column.T*s
    assert sp.expand(z[:2, :2].det()-s.det()) == 0
    assert 4 > 3*1
    assert sp.expand(z.det()**2-s.det()**2*(s[0, 1]-s[1, 0])**2) == 0


if __name__ == "__main__":
    verify()
    print("PASS: nine-word determinant including a regular nilpotent generator and the triple-line commutator cubic")
