"""Check the two determinant models and the universal cubic successor pencil."""
import sympy as sp


def main():
    a, b, c, d = sp.symbols('a b c d')
    coordinates = (a, b, c, d)
    full = sp.Matrix([[a, b], [c, d]])
    trace_free = sp.Matrix([[a, b], [c, -a]])
    assert [sp.hessian(m.det(), coordinates).rank() for m in (full, trace_free)] == [4, 3]
    hankel = sp.Matrix([[a, b], [b, c], [c, d]])
    quadrics = [hankel[[i, j], :].det() for i, j in [(0, 1), (0, 2), (1, 2)]]
    s, t = sp.symbols('s t')
    cubic = dict(zip(coordinates, [s**3, s*s*t, s*t*t, t**3]))
    assert all(sp.expand(q.subs(cubic)) == 0 for q in quadrics)
    u, v, w = sp.symbols('u v w')
    determinant = sp.hessian(sum(x*q for x, q in zip((u, v, w), quadrics)), coordinates).det()
    assert sp.factor(determinant-(u*w-v*v)**2) == 0
    print('PASS: rank-three and rank-four determinants and twisted cubic successor net')


if __name__ == '__main__':
    main()
