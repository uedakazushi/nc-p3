"""Check the matrix covariants, critical derivatives and even-block orientations."""
import sympy as sp

from verify_split_inverse_windows import BITS, cut, reciprocal


def reciprocal_covariants():
    # A full nonsymmetric binary four-tensor checks the adjoint conventions independently of Bell form.
    p = {q: (1+sum((2*i+1)*v for i,v in enumerate(q)))**3+int(len(set(q))==1) for q in BITS}
    eps = sp.Matrix([[0,1],[-1,0]])
    metric = sp.kronecker_product(eps,eps)
    m,n = cut(p,0),cut(p,1)
    gram = metric*m.T*metric*m
    trace,d = gram.trace(),m.det()
    assert 16*n.det() == trace**2-2*(gram**2).trace()+8*d
    f = trace*gram-2*gram**2+2*d*sp.eye(4)
    assert m*gram.inv()*f == 4*n.det()*metric*cut(reciprocal(p,1),0)*metric


def repeated_eigenvalues():
    r,s,z = sp.symbols('r s z',nonzero=True)
    def polynomials(values,d):
        f = sum(values)*z-2*z*z+2*d
        return f, f*f/z
    values = (1,r*s,r,s)
    c = (r-1)*(s-1)
    f,g = polynomials(values,-r*s)
    for v,sign in zip(values,(-1,-1,1,1)):
        assert sp.factor(f.subs(z,v)-sign*c*v)==0
    assert sp.factor(g.diff(z).subs({z:r,s:r})-(r-1)**4)==0
    assert sp.factor(g.diff(z).subs({z:1,s:1/r})+(r-1)*(1/r-1)*(r+1/r-2))==0
    f,g = polynomials((1,1,r,r),r)
    for v in (1,r):
        assert sp.factor(f.subs(z,v)-4*r)==0
        assert sp.factor(g.diff(z).subs(z,v)+16*r)==0
    for values,d in (((1,1,r,r),r),((1,1,-1,-1),-1),((1,1,1,1),1)):
        assert d==values[0]*values[2]


if __name__ == '__main__':
    reciprocal_covariants()
    repeated_eigenvalues()
    print('Inverse-gradient conventions and all repeated-eigenvalue derivatives and determinant signs verified.')
