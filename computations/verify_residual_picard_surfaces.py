"""Mukai lattice and polarization checks for the residual relative Jacobians."""
import sympy as sp


def main():
    t, v, s = sp.symbols('t v s', integer=True)

    def square(r, u):
        return 4*u*u + 16*u*v - 2*r*s

    assert sp.expand(square(-2*t, t) - 4*t*(t+4*v+s)) == 0
    assert sp.expand(square(2*t, 3*t) - 4*t*(9*t+12*v-s)) == 0
    # The two fiber directions in the transformed Mukai lattice.
    pairing = sp.Matrix([[0, 0, 0, -1], [0, 4, 8, 0],
                         [0, 8, 0, 0], [-1, 0, 0, 0]])
    u = sp.Matrix([0, 0, 0, 1])
    root = sp.Matrix([1, 0, 0, 1])
    hyperplane = sp.Matrix([1, 1, 0, 3])

    def reflect(x, vector):
        return x + (x.T*pairing*vector)[0]*vector

    other = reflect(reflect(u, hyperplane), root)
    assert other == sp.Matrix([2, -1, 0, 1])
    assert (u.T*pairing*other)[0] == -2
    assert ((u-other).T*pairing*(u-other))[0] == 4
    print('PASS: both perpendicular lattices, fiber intersection and polarization')


if __name__ == '__main__':
    main()
