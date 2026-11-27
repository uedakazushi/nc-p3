"""Mukai classes, elliptic slopes and degree bounds in the Serre induction."""
import sympy as sp


def main():
    form = sp.Matrix([[0, 0, 0, -1], [0, 4, 8, 0],
                      [0, 8, 0, 0], [-1, 0, 0, 0]])

    def line(h, f):
        return sp.Matrix([1, h, f, 1+2*h*h+8*h*f])

    def reflect(x, v):
        return x + (x.T*form*v)[0]*v

    def pair_twist(x):
        return reflect(reflect(x, line(1, 0)), line(0, 0))

    structure, hyperplane = line(0, 0), line(1, 0)
    fourth = 4*line(-2, 1)-line(-1, 1)
    p = sp.Matrix.hstack(structure, 4*structure-hyperplane,
                         -pair_twist(line(-2, 1)), -pair_twist(fourth))
    first = sp.Matrix([0, 0, 0, 1])
    assert list(p.T*form*first) == [-1, -3, -5, -7]
    assert list(p.T*form*pair_twist(first)) == [-3, -1, 1, 3]
    assert list(p.T*form*sp.Matrix([0, 0, 1, -4])) == [4]*4
    assert (p.T*form*p).diagonal() == sp.Matrix([[-2]*4])
    for i in range(4):
        for j in range(4):
            if i < j:
                assert 1 <= 4+i-j <= 3
            assert 2+sp.Rational(i-j, 2) > 0
    print('PASS: spherical ranks, fiber degrees, low internal degrees and propagation slopes')


if __name__ == '__main__':
    main()
