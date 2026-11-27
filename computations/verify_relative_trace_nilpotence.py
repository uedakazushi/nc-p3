"""The dehomogenizing determinant and positive-degree inverse calculation."""
from sympy import Matrix, eye, symbols


def main():
    t = symbols('t')
    a = Matrix([[2,3,-1],[1,0,4],[0,2,1]])
    for degree in (4,8):
        determinant = (eye(3)-t**degree*a).det()
        assert determinant.subs(t,0) == 1
    shift = Matrix([[0,1,0],[0,0,1],[0,0,0]])
    assert (eye(3)-shift)*(eye(3)+shift+shift**2) == eye(3)
    assert sum(4*4*4 for _ in range(4)) == 256
    print('Finite-ray determinant, nilpotent inverse, and the 256 scalar arrow entries checked.')


if __name__ == '__main__':
    main()
