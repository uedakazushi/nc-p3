"""Check the Euler form and integral Mukai actions for a quartic K3 cut."""
import sympy as sp


def main():
    # (rank, coefficient of h, coefficient of the fiber F, last component).
    pairing = sp.Matrix([[0, 0, 0, -1], [0, 4, 8, 0],
                         [0, 8, 0, 0], [-1, 0, 0, 0]])
    identity = sp.eye(4)

    def line(h, fiber):
        return sp.Matrix([1, h, fiber, 1 + 2*h*h + 8*h*fiber])

    def tensor(h, fiber):
        return sp.Matrix([[1, 0, 0, 0], [h, 1, 0, 0], [fiber, 0, 1, 0],
                          [2*h*h + 8*h*fiber, 4*h + 8*fiber, 8*h, 1]])

    def twist(v):
        return identity + v*v.T*pairing

    objects = sp.Matrix.hstack(line(0, 0), line(1, 0), -line(-2, 1), -line(-1, 1))
    assert -objects.T*pairing*objects == sp.Matrix(
        [[2, 4, 6, 4], [4, 2, 4, 6], [6, 4, 2, 4], [4, 6, 4, 2]])
    first = twist(objects[:, 0])*twist(objects[:, 1])
    product = first*twist(objects[:, 2])*twist(objects[:, 3])
    assert product == first*tensor(-4, 1)*first.inv()*tensor(0, -1)
    degree_four = sp.Matrix([0, 0, 1, -4])
    assert product*degree_four == degree_four
    assert first.inv()*degree_four == sp.Matrix([0, -4, 1, -12])
    print('PASS: Euler form, two-pencil product and degree-four fixed class')


if __name__ == '__main__':
    main()
