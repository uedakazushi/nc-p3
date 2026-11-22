"""Exact quadratic extension and both adjacent-rank conditions on the elliptic chart."""
from itertools import product
from verify_klein_period_four_center import K, paths, commuting_space, artinian_quotient
import sympy as sp
from flint_rational_function import rational_function_field
def numerator(a,b,c,d,z,l):
 return [a+d-b/l-l*c,l*z*a+d/l-z*b-c,a/l+l*z*d-z*b-c,z*(a+d)-l*z*z*b-c/l]
def coeffs(a,b,c,d,z,l):
 D=a*d-b*c
 t=[u/D for u in numerator(a,b,c,d,z,l)]
 f0=(a+d)/(z+1);g0=z*f0
 f1=(t[0]-t[3])/(z*z-1);g1=-t[0]-f1
 e0=h0=(z*b+c)/(z+1);e1=(z*t[1]-t[2])/(z*z-1);h1=(t[2]-e1)/z
 return (f0,f1,g0,g1,e0,e1,h0,h1),(f1*g1,f0*g1+g0*f1+1,f0*g0-a*d),(e1*h1,e0*h1+h0*e1-1,e0*h0-b*c)
def audit(vals, other_sign=False):
 a,b,c,d,z,l=vals
 co,P,Q=coeffs(*vals);f0,f1,g0,g1,e0,e1,h0,h1=co
 p,m,n=P;q,k,o=Q
 assert m*m-4*p*n==k*k-4*q*o
 y0,y1=((m-k)/(2*q),p/q) if other_sign else (-(m+k)/(2*q),-p/q)
 # Each value is represented by its constant and linear coefficient in x.
 unit=a/a;zero=a-a
 def plus(v,w):return tuple(i+j for i,j in zip(v,w))
 def scale(v,k):return tuple(i*k for i in v)
 def mul(v,w):
  v0,v1=v;w0,w1=w;cross=v1*w1
  return (v0*w0-cross*n/p,v0*w1+v1*w0-cross*m/p)
 xv=(zero,unit);yv=(y0,y1);dv=(y0-a*d+b*c,y1+1)
 e=(e0+e1*y0,e1*y1);h=(h0+h1*y0,h1*y1);f=(f0,f1);g=(g0,g1)
 assert plus(mul(e,h),(-b*c-y0,-y1))==(0,0)
 def numer(A,B,C,D,t):
  return list(zip(*(numerator(A[i],B[i],C[i],D[i],z,t) for i in (0,1))))
 N=numer(e,f,g,h,l);N0=numer((a,zero),f,g,(d,zero),unit);N1=numer(e,(b,zero),(c,zero),h,unit)
 for i in range(4):
  det,vec=(yv,N1) if i in (0,3) else (xv,N0)
  lhs,rhs=mul(det,N[i]),mul(dv,vec[i])
  if other_sign and i==0:
   assert lhs!=rhs
   return
  assert lhs==rhs,i
 print('Universal common discriminant and four remaining Schur identities.',flush=True)
def gram_twist(A,B,z,l):
 M=sp.Matrix
 X=M([[1,1/l],[1/l,z]]);Y=M([[1,l],[l*z,1]])
 m0=A.row_join(X).col_join(Y.row_join(Y*A.inv()*X))
 m1=B.row_join(l*X).col_join((Y/l).row_join(Y*B.inv()*X)).T
 def thin(i,j):
  return z if i==3 and j==1 else 1
 def value(i,j,k,n):
  h=i^j
  if h!=k^n:return 0
  if h<2:return (m0,m1)[h][i,k]
  return thin(i,j)*thin(k,n^1)
 w=M(16,16,lambda i,j:value(i//4,i%4,j//4,j%4))
 opp=M(16,16,lambda i,j:value(i//4,j//4,i%4,j%4))
 rev=M(16,16,lambda i,j:value(j//4,i%4,i//4,j%4))
 return w,opp,rev

def four_relations(w,rev):
 out=[]
 flip=[4*j+i for i,j in product(range(4),repeat=2)]
 for p,mat in enumerate((w,rev,w.T,rev.T)):
  rr=sp.Matrix.hstack(*mat.columnspace()).T
  if p%2:rr=rr[:,flip]
  out.append([(i^j,row) for row in rr.tolist() for i,j in [divmod(next(k for k,x in enumerate(row) if x),4)]])
 return out

def rational_example():
    a=sp.Matrix([[1,1],[1,2]])
    b=sp.Matrix([[-sp.Rational(13,14),sp.Rational(1,13)],[-14,-sp.Rational(16,11)]])
    w,opp,rev=gram_twist(a,b,sp.Integer(-2),sp.Integer(-2))
    assert w.rank()==rev.rank()==6
    assert opp.det()==-sp.Rational(29506135761,19712)
    vals=tuple(map(sp.Integer,(1,1,1,2,-2,-2)))
    co,p,q=coeffs(*vals)
    x=a[0,0]*a[1,1]-b[0,1]*b[1,0]
    y=b[0,0]*b[1,1]-a[0,1]*a[1,0]
    assert 2*p[0]*x+p[1]==-2*q[0]*y-q[1]==-sp.Rational(27,2)
    assert p[0]*x*x+p[1]*x+p[2]==q[0]*y*y+q[1]*y+q[2]==0
    f0,f1,g0,g1,e0,e1,h0,h1=co
    assert b==sp.Matrix([[e0+e1*y,f0+f1*x],[g0+g1*x,h0+h1*y]])
    return w,opp,rev


if __name__=='__main__':
    names=('a','b','c','d','z','l')
    field=rational_function_field(names)
    audit([field.from_sympy(v) for v in sp.symbols(names)])
    audit(tuple(map(sp.Integer,(5,2,3,7,-2,2))),other_sign=True)
    l=sp.Symbol('l')
    co,p,q=coeffs(*tuple(map(sp.Integer,(5,2,3,7,-2))),l)
    quartic=sp.cancel(l*l*(p[1]**2-4*p[0]*p[2]))
    assert sp.Poly(quartic,l).degree()==4 and sp.discriminant(quartic,l)!=0
    w,opp,rev=rational_example()
    relations=four_relations(w,rev)
    basis,image=paths(relations,5)
    assert all([sum(len(basis[p,n,h]) for h in K) for n in range(6)]==[1,4,10,20,35,56] for p in K)
    assert [commuting_space(basis,image,h,0).cols for h in K]==[3,2,0,0]
    artinian_quotient(relations,basis,image)
    print('Smooth elliptic parameter curve, rank sixteen and rational central-parameter certificate.')
