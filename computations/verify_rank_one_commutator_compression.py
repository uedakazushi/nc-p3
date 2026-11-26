"""Rank-one commutator images, orbit forms, and the incident flag equations."""
import sympy as sp


def comm(x, y):
    return x*y-y*x


def check_image(p, positions):
    units = []
    for i in range(3):
        for j in range(3):
            e = sp.zeros(3)
            e[i, j] = 1
            units.append(e)
    images = [comm(p, x) for x in units]
    assert all(not any(z[i, j] for i, j in positions) and z.trace() == 0 for z in images)
    mat = sp.Matrix.hstack(*(z.reshape(9, 1) for z in images))
    _, pivots = mat.rref()
    assert len(pivots) == 4
    forms = sp.Matrix(4, 4, lambda i, j: (p*comm(units[pivots[i]], units[pivots[j]])).trace())
    assert forms == -forms.T and forms.det() != 0


def verify():
    incident = sp.zeros(3)
    incident[0, 2] = 1
    check_image(incident, [(1, 0), (1, 1), (2, 0), (2, 1)])
    check_image(sp.diag(0, 0, 1), [(0, 0), (0, 1), (1, 0), (1, 1)])
    x = sp.Matrix(3, 3, sp.symbols("x:9"))
    y = sp.Matrix(3, 3, sp.symbols("y:9"))
    def wedge(i, j):
        return x[i]*y[j]-x[j]*y[i]
    bracket = comm(x, y)
    expressions = [wedge(3, 0)-wedge(3, 4)+wedge(5, 6),
                   wedge(6, 0)-wedge(6, 8)+wedge(7, 3),
                   wedge(6, 1)+wedge(7, 4)-wedge(7, 8),
                   wedge(3, 1)+wedge(5, 7)]
    j = sp.Symbol("j", nonzero=True)
    a, h, d, e = sp.symbols("a h d e")
    metric = sp.Matrix([[0, 0, 1], [0, j, 0], [1, 0, 0]])
    compressed = sp.Matrix([[a, h, e], [0, 0, d], [0, 0, -a]])
    constraints = metric*compressed+(metric*compressed).T
    assert constraints.subs({e: 0, h: -j*d}) == sp.zeros(3)
    assert sp.linear_eq_to_matrix(list(constraints), (a, h, d, e))[0].rank() == 2
    for pos, expression in zip([(1, 0), (2, 0), (2, 1), (1, 1)], expressions):
        assert sp.expand(bracket[pos]-expression) == 0


if __name__ == "__main__":
    verify()
    print("PASS: both four-dimensional compression images, orbit symplectic forms, and flag equations")
