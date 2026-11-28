"""Identities for the smooth quadric reshaping and the cone's binary-cubic pencil."""

import sympy as sp


def main():
    a,b,c,d,e,f,g,h,i = sp.symbols('a b c d e f g h i')
    p = sp.Matrix([[a,b,d,e],[b,c,e,f],[d,e,g,h],[e,f,h,i]])
    flat = sp.Matrix([[a,b,c],[d,e,f],[g,h,i]])
    adj = p.adjugate()
    assert sp.expand(adj[0,3]-adj[1,2]-flat.det()) == 0
    entries = sp.symbols('u:6')
    a,b,c,d,e,f = entries
    block = sp.Matrix([[a,b,c],[b,d,e],[c,e,f]])
    column = sp.Matrix(sp.symbols('v:3'))
    scalar = sp.symbols('s')
    schur = scalar*block-column*column.T
    def klein(m):
        return m[0,0]*m[2,2]-m[0,1]*m[1,2]+m[0,2]*(m[1,1]-m[0,2])
    identity = sp.Poly(klein(schur),scalar)
    assert identity.nth(0) == 0 and sp.expand(identity.nth(2)-klein(block)) == 0
    coefficients = sp.hessian(identity.nth(1), list(column))
    polar = sp.Matrix([coefficients[i,j] for i in range(3) for j in range(i,3)]).jacobian(entries)
    assert polar.det() != 0
    assert sp.expand(klein(schur)-scalar*(block.row_join(column).col_join(column.T.row_join(sp.Matrix([[scalar]]))).adjugate()[1,1]
                       -block.row_join(column).col_join(column.T.row_join(sp.Matrix([[scalar]]))).adjugate()[0,2])) == 0
    print('Verified the moment determinant identity and the nonsingular six-coordinate cone polarization.')


if __name__ == '__main__':
    main()
