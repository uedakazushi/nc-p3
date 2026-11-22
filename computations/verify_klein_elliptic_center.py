"""Exact quartic centers on the smooth discriminant curve over Q(lambda)."""
from flint import fmpq_poly as Poly,fmpq
from verify_klein_elliptic_completion import coeffs, rational_example, four_relations
from verify_klein_center_lift import field_paths,central_equations,echelon,right_inverse_product,transpose
import sympy as sp
from itertools import product
def curve_field(quartic):
 one=Poly([1]);zero=Poly([])
 class E:
  __slots__=('u','v','d')
  def __init__(self,u,v=zero,d=one):
   if not d:raise ZeroDivisionError
   if d!=one:
    if d.degree()>0:
     g=d.gcd(u).gcd(v)
     if g.degree()>0:u,v,d=u//g,v//g,d//g
    lc=d.leading_coefficient();u,v,d=u/lc,v/lc,d/lc
   self.u,self.v,self.d=u,v,d
  @staticmethod
  def co(x):return x if isinstance(x,E) else E(Poly([x]))
  def __bool__(self):return bool(self.u) or bool(self.v)
  def __eq__(self,o):
   o=self.co(o);return self.u==o.u and self.v==o.v and self.d==o.d
  def __add__(self,o):
   o=self.co(o)
   if not o:return self
   if not self:return o
   if o.d==self.d:return E(self.u+o.u,self.v+o.v,self.d)
   return E(self.u*o.d+o.u*self.d,self.v*o.d+o.v*self.d,self.d*o.d)
  __radd__=__add__
  def __neg__(self):return E(-self.u,-self.v,self.d)
  def __sub__(self,o):return self+-self.co(o)
  def __rsub__(self,o):return self.co(o)+-self
  def __mul__(self,o):
   o=self.co(o)
   return E(self.u*o.u+quartic*self.v*o.v,self.u*o.v+self.v*o.u,self.d*o.d)
  __rmul__=__mul__
  def __truediv__(self,o):
   o=self.co(o)
   if not o.v:return E(self.u*o.d,self.v*o.d,self.d*o.u)
   return E((self.u*o.u-quartic*self.v*o.v)*o.d,(self.v*o.u-self.u*o.v)*o.d,self.d*(o.u*o.u-quartic*o.v*o.v))
  def __rtruediv__(self,o):return self.co(o)/self
 return E
def rels(a,b,c,d,z,l,root):
 co,P,Q=coeffs(a,b,c,d,z,l)
 p,m,n=P;q,k,o=Q
 assert root*root==m*m-4*p*n==k*k-4*q*o
 x=(root-m)/(2*p);y=(-root-k)/(2*q)
 f0,f1,g0,g1,e0,e1,h0,h1=co;e,h=e0+e1*y,h0+h1*y;f,g=f0+f1*x,g0+g1*x
 A=[[a,b],[c,d]];B=[[e,f],[g,h]];C=[[a,f],[g,d]];D=[[e,b],[c,h]]
 X=[[1,1/l],[1/l,z]];Y=[[1,l],[l*z,1]]
 XX=[[1,1],[1,z]];YY=[[1,1],[z,1]]
 def times(t,M):return [[t*u for u in row] for row in M]
 rs=[]
 for cut in range(4):
  rev=cut%2;dual=cut//2
  if rev:S,T,V,W,rr=C,D,XX,YY,l/l
  else:S,T,V,W,rr=A,B,X,Y,l
  graphs=(right_inverse_product(W,S),right_inverse_product(transpose(times(rr,V)),transpose(T))) if not dual else (right_inverse_product(transpose(V),transpose(S)),right_inverse_product(times(1/rr,W),T))
  rows=[]
  for h,graph in enumerate(graphs):
   for i in range(2):
    coeff=[int(i==j) for j in range(2)]+[graph[j][i] for j in range(2)]
    rows.append((h,{(j,j^h):val for j,val in enumerate(coeff) if val}))
  for h in (2,3):
   coeff=[z if i==3 and (i^h^dual)==1 else 1 for i in range(4)]
   if rev and h==3:coeff=[u*(1/l if dual else l) if i<2 else u for i,u in enumerate(coeff)]
   rows.append((h,{(i,i^h):u for i,u in enumerate(coeff)}))
  rs.append([(h,{(j,i):t for (i,j),t in row.items()}) for h,row in rows] if rev else rows)
 return rs
def check_graphs():
 w,_,n=rational_example()
 vals=tuple(map(sp.Rational,(1,1,1,2,-2,-2,-sp.Rational(27,2))))
 for original,graph in zip(four_relations(w,n),rels(*vals)):
  a=sp.Matrix([r for _,r in original])
  b=sp.Matrix([[r.get(ij,0) for ij in product(range(4),repeat=2)] for _,r in graph])
  assert a.col_join(b).rank()==b.rank()==6

if __name__=='__main__':
 check_graphs()
 coeff=[7756,8856,-22355,-17712,31024];quartic=Poly(coeff)/fmpq(7569)
 E=curve_field(quartic);l=E(Poly([0,1]));r=E(Poly([]),Poly([1]))/l
 params=list(map(E.co,[5,2,3,7,-2]))+[l,r]
 R=rels(*params);bs,im=field_paths(R,E.co(1));assert all([sum(len(bs[p,n,h]) for h in range(4)) for n in range(6)]==[1,4,10,20,35,56] for p in range(4))
 for h in (0,1):
  rows,dim=central_equations(bs,im,h);dimension=dim-len(echelon(rows))
  assert dimension==[3,2][h]
  print('Genus-one function field, character',h,'quartic center dimension',dimension,flush=True)
