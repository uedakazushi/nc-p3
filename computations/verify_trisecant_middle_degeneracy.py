"""Rank-two Sylvester middle, basepoint-free quintics, and the twenty-dimensional cubic quotient."""

import sympy as sp


def main():
    columns = []
    for i in range(5):
        for j in range(5):
            block = sp.zeros(6,4)
            if j < 4:
                block[i+1,j] += 1
            if j:
                block[i,j-1] -= 1
            columns.append(block.reshape(24,1))
    operator = sp.Matrix.hstack(*columns)
    assert operator.rank() == 24
    desired = sp.zeros(6,4)
    # The two interior coefficient functionals span a line avoiding all pure cubics.
    desired[0,1],desired[5,2] = 1,1
    vector, parameters = operator.gauss_jordan_solve(desired.reshape(24,1))
    transform = vector.subs({parameters[0]:3}).reshape(5,5)
    assert transform.det() and operator*transform.reshape(25,1) == desired.reshape(24,1)
    assert desired.rank() == 2
    # Chi of the internal degree-three Koszul complex on P^3 of the middle.
    t = sp.symbols('t')
    chi = lambda k: sp.prod(t+k+i for i in (1,2,3))/6
    # Source O(-1)^4 dualizes to O(1)^4; the target is trivial rank six.
    length = 20*chi(3)-6*10*chi(2)+15*4*chi(1)-20*chi(0)
    assert sp.expand(length) == 20
    print('Verified the commutator quotient, explicit invertible rank-two example, and Koszul length twenty.')


if __name__ == '__main__':
    main()
