"""Orthogonal substitutions, curve parameters and rotations on the scalar Jordan spectrum."""
import sympy as sp

from verify_split_inverse_windows import BITS, cut, positions
from verify_split_two_eigenvalue_jordan import diagonal_invariants

T, U, Z = sp.symbols('t u z')


def chain(m, t):
    return sp.Matrix(m,m,lambda i,j:sp.binomial(sp.Rational(m-2,2)-j,i-j)*t**(i-j)
                     if i>=j else 0)


def verify_chains():
    for m in range(1,5):
        n = sp.zeros(m)
        for j in range(m-1):
            n[j+1,j] = 1
        identity = sp.eye(m)
        gram = identity[:,::-1]
        change = chain(m,T)
        assert (change.T*gram*change-gram).applyfunc(sp.expand).is_zero_matrix
        assert (change*chain(m,U)-chain(m,T+U)).applyfunc(sp.expand).is_zero_matrix
        assert (change*n-n*(identity+T*n).inv()*change).applyfunc(sp.expand).is_zero_matrix
        k = n/2-n**2/8+n**3/16
        half = chain(m,sp.Rational(1,2))
        assert half*(identity-k)*half.inv() == (identity+k).inv()
        assert (identity+k).inv()*(identity-n**2/2) == identity-n/2-n**2/8-n**3/16
    change = sp.diag(chain(3,T),1)
    f = sp.zeros(4)
    f[3,0],f[2,3] = 1,-1
    reflection = sp.diag(1,-1,1,1)
    assert reflection.det() == -1
    assert change*f == f*change and reflection*f*reflection == f
    assert reflection*change*reflection == sp.diag(chain(3,-T),1)


def verify_rotation():
    j = sp.Matrix([[0,-1],[1,0]])
    metric = sp.kronecker_product(j,j)
    switch = sp.eye(4)[:,[0,2,1,3]]
    e = sp.Matrix([[0,1],[0,0]])
    z = sp.diag(1,-1)
    nil = sp.kronecker_product(e,z)+sp.kronecker_product(z,e)
    mat = sp.eye(4)+nil
    raw = mat*metric
    w = {positions(a,b,0):raw[a,b] for a in range(4) for b in range(4)}
    rotated = {q[1:]+q[:1]:v for q,v in w.items()}
    assert cut(rotated,0)*metric == -switch*(sp.eye(4)-nil)*switch
    inverse = mat.inv()*metric
    terms = {positions(a,b,0):inverse[a,b] for a in range(4) for b in range(4)}
    swapped = {q[::-1]:v for q,v in terms.items()}
    assert cut(swapped,0)*metric == switch*mat.inv()*switch
    assert switch.det() == -1


if __name__ == '__main__':
    verify_chains()
    verify_rotation()
    diagonal_invariants(sp.Integer(1), sp.Symbol('a',nonzero=True), None)
    print('All chain substitutions, the one-shift, exchange on the 3+1 parameter, and the toral row verified.')
