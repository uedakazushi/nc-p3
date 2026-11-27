"""Apolar quartics, uniform negative point ideals and period-one tests on the smooth chart."""
from itertools import combinations
import sympy as s
from verify_klein_elliptic_completion import gram_twist
from verify_negative_component_point_graph import K,D,a,b,c,d,z,x,evaluation,quadrics


def apolar():
    A = s.Matrix([[a,-b],[-c,d]])
    W,_,_ = gram_twist(A,A,z,s.S.One)
    blocks = [s.Matrix(4,4,lambda i,j:W[4*i+j,4*(i^h)+(j^h)]) for h in K]
    factors = (a-d,b*z-c,a+b,a+c,b+d,c+d,a+b*z,a*z+c,c+d*z,b*z+d)
    for h,m in enumerate(blocks):
        assert m.rank()==2
        if h<2:
            n = (m.T if h else m)[:,[0,2]]
            minors = [s.factor(n[list(J),:].det()) for J in combinations(K,2)]
            product = s.factor(s.prod(minors))
            assert s.factor(product**2*D**10/((z-1)**2*s.prod(factors)**2))==1
        else:
            assert any(s.factor(m.extract(J,L).det()) in (z-1,1-z)
                       for J in combinations(K,2) for L in combinations(K,2))


def point_ideals():
    E = evaluation()
    def minor(i):
        return E.extract(tuple(map(int,i)),K).det(method='domain-ge')
    X,Y,Z,T = x
    P = D*(X**2-Y**2)-(a+b+c+d)*Z**2+(a*z+b*z*z+c+d*z)*T**2
    Q = a*X**2-d*Y**2-Z**2+z*T**2
    tests = {'0145':X*Y*P+(z-1)*Z*T*Q,'2345':-X*Y*P+(z-1)*Z*T*Q,
             '0124':-(a+b*z)*X*T*Q+Y*Z*(P-(a+b)*Q),
             '0125':-(a+b)*X*Z*Q+Y*T*(P-(a+b*z)*Q),
             '0134':X*T*(z*P-(c+d*z)*Q)-(c+d)*Y*Z*Q,
             '0135':X*Z*(P-(c+d)*Q)-(c+d*z)*Y*T*Q,
             '0234':(a+c)*X*Z*Q+Y*T*(-z*P+(a*z+c)*Q),
             '0235':(a*z+c)*X*T*Q+Y*Z*(-P+(a+c)*Q),
             '1234':X*Z*(-P+(b+d)*Q)+(b*z+d)*Y*T*Q,
             '1235':X*T*(-P+(b*z+d)*Q)+(b+d)*Y*Z*Q,
             '1245':-Q*(b*(X**2-Y**2)+Z**2-T**2)}
    assert all(s.expand(minor(I)-f)==0 for I,f in tests.items())
    eq = quadrics(False)
    M = s.Matrix([P,Q]).jacobian((X,Y))*s.diag(1/(2*X),1/(2*Y))
    C = s.Matrix(eq).jacobian((X,Y))*s.diag(1/(2*X),1/(2*Y))
    assert all(s.factor(f)==0 for f in s.Matrix(eq)-C*M.inv()*s.Matrix([P,Q]))
    assert s.factor((a+b)*(c+d*z)-(a+b*z)*(c+d))==(z-1)*D


def group_on_fiber():
    p,q,r,t = s.symbols('p q r t')
    B = s.Matrix([[2*z,z+1],[-z*(z+1),-2*z]])
    J = s.Matrix([[0,1],[z,0]])
    assert B*J+J*B==s.zeros(2)
    assert (B**2).applyfunc(s.factor)==-z*(z-1)**2*s.eye(2)
    for r0,sign in ((-p,-1),(-z*p,1)):
        v = s.Matrix([p,r0])
        assert ((z-1)*J*v-sign*B*v).applyfunc(s.factor)==s.zeros(2,1)
    coords = {a:(p+q)/2,d:(p-q)/2,b:(r+t)/(2*z),c:(r-t)/2}
    assert s.factor(4*D.xreplace(coords)-(p*p-q*q-(r*r-t*t)/z))==0


if __name__=='__main__':
    apolar()
    point_ideals()
    group_on_fiber()
    print('PASS: smoothness divisor, uniform projected point schemes and both period-one planes')
