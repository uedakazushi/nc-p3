"""Seven-dimensional involution locus from four independent periodic point tracks."""
from functools import lru_cache
from itertools import product

import numpy as np
import sympy as sp
from flint import nmod_mat

from flint_rational_function import rational_function_field
from verify_exchange_mixed_configurations import modular_points, partner
from verify_klein_elliptic_completion import rational_example
K=range(4);allow=[w for w in product(K,repeat=4) if all(w[i]!=w[(i+1)%4] for i in K)]
def rows(A,p):return nmod_mat([[int(x)%p for x in y] for y in A],p)
def kernel(A):
 B,k=A.nullspace();return np.array(B.tolist(),dtype=np.int64)[:,:k].T
def jac(T,p,orbits):
 jac=[]
 for off in (0,1):
  perm=tuple((off+j)%4 for j in K)
  M=T.transpose(perm).reshape(16,16);assert rows(M,p).rank()==6
  L,R=kernel(rows(M.T,p)),kernel(rows(M,p))
  J=np.array([sum(np.outer(L[:,4*w[perm[0]]+w[perm[1]]],R[:,4*w[perm[2]]+w[perm[3]]]) for w in orbit).ravel()%p for orbit in orbits]).T
  jac.extend(J)
 # the 2-by-2 opposite determinants for cross-pair labels
 extra=[]
 for off in (0,1):
  Q=T.transpose(tuple((off+i)%4 for i in K))
  for a,c in product(K,repeat=2):
   if a//2==c//2:continue
   i,j=[k for k in K if k!=a and k!=c]
   base=Q[a,:,c,:][np.ix_([i,j],[i,j])]
   assert (base[0,0]*base[1,1]-base[0,1]*base[1,0])%p==0
   grad={}
   for x,y in product((0,1),repeat=2):
    q=(a,(i,j)[x],c,(i,j)[y]);q=q[-off:]+q[:-off] if off else q
    grad[q]=(-1)**(x+y)*base[1-x,1-y]
   extra.append([sum(grad.get(w,0) for w in orb)%p for orb in orbits])
 r=rows(jac,p).rank()
 return len(orbits)-r,rows(jac+extra,p).rank()-r
def sample():
 w,_,n=rational_example();p=2657;xx=[modular_points(w,p)]
 for mat,t in [(w,False),(n,True),(w.T,False)]:xx.append(partner(xx[-1],mat,t,p))
 Q=np.array([int(v.p)*pow(int(v.q),-1,p)%p for v in w]).reshape([4]*4)
 frame=[2,17,5,14]
 X=np.array([[s for s in (np.array(x)[frame])] for x in xx])
 # normalize the pairs by one common first coordinate (all equal to 1)
 for t in range(4):Q=np.moveaxis(np.tensordot(X[t],Q,axes=(1,t))%p,0,t)
 assert all(Q[w]==0 for w in product(K,repeat=4) if w not in allow)
 assert all(Q[w]!=0 for w in allow)
 return Q,p

dot=[(0,1),(2,3),(0,2),(0,3),(2,0),(2,1)]
def up(q):return tuple(i^1 for i in q) if q[0]%2 else q
def coeff(w):return sp.Symbol('w'+''.join(map(str,up(w)))) if w in allow else 0
def matrix(shift,sign):
 def entry(i,j,k,l):
  n=(i,j,k,l)
  return coeff(n[-shift:]+n[:-shift] if shift else n)
 return sp.Matrix(6,6,lambda i,j:entry(*dot[i],*dot[j])+sign*entry(*dot[i],*(x^1 for x in dot[j])))

eqs=set()
for off in range(2):
 for a,c in product(K,repeat=2):
  if a//2==c//2:continue
  i,j=[k for k in K if k not in (a,c)]
  val=lambda x,y: coeff((a,x,c,y) if off==0 else (y,a,x,c))
  eq=val(i,i)*val(j,j)-val(i,j)*val(j,i)
  eqs.add(eq if str(eq).startswith('w') else -eq)

corr={}

for eq in eqs:
 elim=next(q for q in eq.free_symbols if str(q) in ('w2023','w2123','w0323','w0232','w2320','w2321','w2313','w2303'))
 corr[elim]=sp.solve(eq,elim)[0]
renames={'0101':'P','2323':'Q','0123':'A','0132':'B','2301':'C','2310':'D','0231':'E','0321':'F','2013':'G','2103':'H',
         '0102':'i','0103':'j','0121':'k','0131':'l','0201':'m','0301':'n','2010':'o','2101':'r'}
rename={sp.Symbol('w'+w):sp.Symbol(z) for w,z in renames.items()}
def mat(p,sign):return matrix(p,sign).subs(corr).xreplace(rename)

x,xx,u,uu,B,F,k=sp.symbols('x X u U B F k')
corevals=[1,1,1,B,(1-x)*(1-xx),x*xx/B,1,F,u*uu,(1-u)*(1-uu)/F,
          B*u,(1-u)/k,k,1,x,F*(1-x)/k,xx*uu/B,k*(1-xx)*(1-uu)/F]
sub=dict(zip(map(sp.Symbol,renames.values()),corevals))
def norm(M):return M.applyfunc(sp.factor)
def matrices():
 DM=[]
 for p in range(2):
  neg=norm(mat(p,-1).subs(sub));pos=norm(mat(p,1).subs(sub))
  C=norm(neg[2:,:2]*neg[:2,:2].adjugate()*neg[:2,2:]/neg[:2,:2].det())
  D=norm(neg[2:,2:]-C)
  assert D[:2,2:]==sp.zeros(2) and D[2:,:2]==sp.zeros(2)
  DM.append((pos,D))
 unknown=sorted({z for _,D in DM for v in D for z in v.free_symbols if str(z).startswith('w')},key=str)
 L,rhs=sp.linear_eq_to_matrix([D[i,j] for _,D in DM for i,j in product(range(4),repeat=2) if i//2==j//2],unknown)
 assert L.shape==(16,16) and L.rank()==12
 obstructions=[sp.factor((z.T*rhs)[0]) for z in L.T.nullspace()]
 piv=L.rref()[1]
 sel=L[:,piv].T.rref()[1]
 zz={unknown[i]:unknown[i] for i in range(len(unknown)) if i not in piv}
 val=L.extract(sel,piv).inv()*(rhs[list(sel),:]-L.extract(sel,[i for i in range(len(unknown)) if i not in piv])*sp.Matrix(list(zz)))
 zz.update({unknown[i]:sp.factor(v) for i,v in zip(piv,val)})
 return [norm(m.subs(zz)) for m,_ in DM],zz,obstructions

field=rational_function_field(('x','X','u','B','F'))
zero,one=field.zero,field.one
env={s:(field.from_sympy(s),zero) for s in (x,xx,u,B,F)}
delta=one
def mul(a,b):return a[0]*b[0]+delta*a[1]*b[1],a[1]*b[0]+a[0]*b[1]
def power(a,n):
 if n<0:
  d=a[0]**2-delta*a[1]**2
  return power((a[0]/d,-a[1]/d),-n)
 if n==0:return one,zero
 if n==1:return a
 y=power(a,n//2);v=mul(y,y)
 return mul(v,a) if n%2 else v
@lru_cache(None)
def ev(e):
 if e in env:return env[e]
 if e.is_Number:return field.from_sympy(e),zero
 if e.is_Add:
  y=[ev(v) for v in e.args];return sum((v[0] for v in y),zero),sum((v[1] for v in y),zero)
 if e.is_Mul:
  p=(one,zero)
  for v in e.args:p=mul(p,ev(v))
  return p
 if e.is_Pow:return power(ev(e.base),int(e.exp))
 raise ValueError(e)
def setup(w):
 global delta

 P=sp.Poly(sp.cancel(sp.fraction(w[0])[0]/k),uu)
 assert P.degree()==1
 env[uu]=(-ev(P.nth(0))[0]/ev(P.nth(1))[0],zero)
 L=sp.Poly(sp.fraction(w[1])[0],k)
 assert L.degree()==2 and L.nth(1)==0
 delta=-ev(L.nth(0))[0]/ev(L.nth(2))[0]
 env[k]=zero,one
 for q in w:assert ev(q)==(zero,zero)
 return P,L
def verify(mats,Z):
 vs={v:0 for v,z in Z.items() if v==z}
 S=[norm(m[2:,2:]-norm(m[2:,:2]*m[:2,:2].adjugate()*m[:2,2:]/m[:2,:2].det())) for m in mats]
 for s in S:assert s[:2,2:]==s[2:,:2]==sp.zeros(2)
 for i,j in ((0,2),(2,0)):
  polys=[];lin=[]
  for block in (S[0][i:i+2,i:i+2],S[1][j:j+2,j:j+2]):
   c=[];vars=[]
   for v in block:
    q=sp.factor(v.subs(vs));c.append(ev(q));vars.append(sp.factor(v-q))
   assert vars[0]==vars[3] and vars[1]==vars[2]
   lin.append(vars)
   av=lambda a,b:(a[0]+b[0],a[1]+b[1])
   ad=mul(c[0],c[3]);bc=mul(c[1],c[2])
   polys.append((av(c[0],c[3]),av(c[1],c[2]),(ad[0]-bc[0],ad[1]-bc[1])))
  assert lin[0]==lin[1] and polys[0]==polys[1]


def residue(v,p):
    n,d=sp.fraction(v)
    return int(n)*pow(int(d),-1,p)%p


def audit_basepoint(P0,L0,mats,Z):
    T,p=sample()
    fixed=[(w,tuple(a^1 for a in w)) for w in allow if w[0]%2==0]
    assert jac(T,p,fixed)==(12,0)
    assert jac(T,p,[(w,) for w in allow])==(28,2)
    C={name:int(T[tuple(map(int,w))]) for w,name in renames.items()}
    div=lambda v,d:v*pow(d,-1,p)%p
    raw=[div(C[t]*C[r],C[z]*C['P']) for t,r,z in [('l','m','E'),('i','o','G'),('i','l','B'),('m','o','D')]]
    spec=dict(zip((x,xx,u,uu,B,F,k),raw+[div(C['B'],C['A']),div(C['F'],C['E']),div(C['k'],C['l'])]))
    assert raw[0]!=raw[1] and raw[2]!=raw[3]
    assert residue(P0.nth(1).subs(spec),p) and residue(L0.nth(2).subs(spec),p)
    # The following diagonal rescaling fixes P,Q,A,E,l to one.
    t2=div(C['P'],C['l']);t3=div(C['P'],C['A']*t2);t1=div(C['P'],C['E']*t2)
    t0=div(C['P'],C['Q']*t1*t2*t3);scales=[t0,t1,t2,t3]
    def entry(w):
        scale=sp.prod(scales[i] for i in K if w[i]>=2)
        return div(int(scale)*int(T[w]),C['P'])
    actual={coeff(w):entry(w) for w in allow}
    for w,z0 in renames.items():
        z=sp.Symbol(z0)
        assert residue(sp.sympify(sub[z]).subs(spec),p)==actual[sp.Symbol('w'+w)]
    for w,f in corr.items():assert residue(f.subs(actual),p)==actual[w]
    for w,f in Z.items():
        assert residue(f.subs(spec).subs(actual),p)==actual[w]
    for off in (0,1):
        for sign,required in ((1,4),(-1,2)):
            mat0=matrix(off,sign).subs(actual)
            assert rows(mat0.tolist(),p).rank()==required
            assert int(mat0[:2,:2].det())%p
        assert norm(mats[off].subs(spec).subs(actual)).applyfunc(lambda z:residue(z,p))==matrix(off,1).subs(actual).applyfunc(lambda z:int(z)%p)


if __name__=='__main__':
    mats,Z,obstructions=matrices()
    P0,L0=setup(obstructions)
    verify(mats,Z)
    audit_basepoint(P0,L0,mats,Z)
    print('PASS: four opposite-cut identities, two conics, full incidence dimension and the marked full-rank specialization')
