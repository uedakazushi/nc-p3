"""Verify the diagonal divisor, the two-gon smoothing and the binary pencil on its normalization."""
import sympy as sp


def main():
    u = sp.Matrix([1,1,1,1])
    v = sp.Matrix([1,-1,0,0])
    w = sp.Matrix([0,0,1,-1])
    ns = 4*sp.eye(4)
    for i in (0,1):
        for j in (2,3):
            delta = sp.eye(4)[:,i]-sp.eye(4)[:,j]
            ns -= delta*delta.T
    assert ns == u*u.T + v*v.T + w*w.T
    assert ns*sp.Matrix([1,1,-1,-1]) == sp.zeros(4,1)
    assert ns[:3,:3].det() == 4
    signs = [(int(b==0),int((b^1)==0)) for b in range(4)]
    assert sorted(signs) == [(0,0),(0,0),(0,1),(1,0)]

    x,y,t = sp.symbols("x y t")
    equation = y**2 - x*(x-1)*(x-t)
    scale = sp.symbols("scale")
    tangent = sp.expand(equation.subs({z:scale*z for z in (x,y,t)})).coeff(scale,2)
    assert sp.expand(equation-tangent) == -x**3+t*x**2
    assert sp.hessian(tangent,(x,y,t)).det() != 0
    assert tangent.subs({x:0,y:0,t:1}) == 0
    assert sp.factor(tangent.subs(t,0),extension=sp.I) == (x-sp.I*y)*(x+sp.I*y)

    T = [2,3,5,7,11,13,17,19,23,29,31]
    Y = [3,1,4,1,5,9,2,6,5]
    G = sp.Matrix(6,6,lambda i,j:T[i+j])
    H = sp.Matrix(6,4,lambda i,j:Y[i+j])
    L = G.inv()*H
    ideal = sp.Matrix(5,7,lambda i,j:T[i+j]).nullspace()
    rows = sp.Matrix(2,3,lambda i,j:sum(ideal[i][k]*Y[j+k] for k in range(7)))
    ell = rows.row(0).cross(rows.row(1)).T
    polar = sp.Matrix(7,5,lambda i,j:T[i+j])
    rhs = sp.Matrix(7,3,lambda i,j:Y[i+j])*ell
    c = polar[:5,:].inv()*rhs[:5,:]
    def polynomial(coeffs):
        return sum(a*x**j for j,a in enumerate(coeffs))
    ell,c = polynomial(ell),polynomial(c)
    assert sp.discriminant(ell,x) and sp.resultant(ell,c,x)
    def column(f):
        return sp.Matrix([sp.expand(f).coeff(x,i) for i in range(8)])
    domain = sp.Matrix.hstack(*(column(ell*x**i) for i in range(6)),*(column(c*x**j) for j in range(4)))
    target = sp.eye(6).row_join(L)
    columns = domain.rref()[1]
    projection = target[:,columns]*domain[:,columns].inv()
    assert projection*domain == target
    pencil = [polynomial(q) for q in projection.nullspace()]
    assert len(pencil) == 2 and sp.gcd(*pencil) == 1
    assert max(sp.degree(q,x) for q in pencil) == 7
    print("Four-diagonal class and parities, Legendre two-gon, and the basepoint-free degree-seven pencil verified.")


if __name__ == '__main__':
    main()
