"""Klein equation for a ternary tangent net and symmetric triple-cone examples."""

import sympy as sp


def bezout(a, b, z, w):
    numerator = (a*b.subs(z, w)-b*a.subs(z, w)).expand()
    quotient, remainder = sp.div(numerator, z-w, z)
    assert not remainder
    poly = sp.Poly(quotient, z, w)
    return sp.Matrix(3, 3, lambda i,j: poly.coeff_monomial(z**i*w**j))


def main():
    a,b,c,d,e,f = sp.symbols('a b c d e f')
    matrix = sp.Matrix([[a,b,c],[b,c+d,e],[c,e,f]])
    adj = matrix.adjugate()
    assert sp.expand(adj[1,1]-adj[0,2]-(a*f-b*e+c*d)) == 0
    z,w,x,y,u = sp.symbols('z w x y u')
    fixed = (z-1)*(z-2)*(z+1)
    moving = x*z**2+y*z+u
    assert sp.expand(bezout(fixed,moving,z,w).det()+sp.resultant(fixed,moving,z)) == 0
    basis = [sp.Integer(1),z**2,z**3-z]
    net = x*bezout(basis[1],basis[2],z,w)+y*bezout(basis[2],basis[0],z,w)+u*bezout(basis[0],basis[1],z,w)
    cubic = net.det().expand()
    assert len(sp.factor_list(cubic)[1]) == 1 and sp.Poly(cubic,x,y,u).total_degree() == 3
    assert cubic.subs({x:1,y:z**2,u:z**3-z}).expand() == 0
    assert all(sp.diff(cubic,v).subs({x:1,y:1,u:0}) == 0 for v in (x,y,u))
    assert net.subs({x:1,y:1,u:0}).rank() == 1
    print('Verified the Klein quadric identity, the three-line cone, and the integral nodal cubic cone.')


if __name__ == '__main__':
    main()
