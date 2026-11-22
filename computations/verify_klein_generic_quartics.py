"""Certify the five central quartics over four independent characteristic-zero parameters."""
import sympy as sp
from flint_rational_function import rational_function_field
from verify_klein_center_lift import field_relations, field_paths, central_equations


def cleared(row):
 if not row:return {}
 els=list(row.values());unit=els[0].d**0
 den=unit
 for a in els:
  if a.d!=unit and a.d!=den:den=(den//den.gcd(a.d))*a.d
 return primitive({k:a.n*(den//a.d) for k,a in row.items()})
def primitive(r):
 if not r:return r
 ordered=sorted(r.values(),key=len)
 g=ordered[0]
 for a in ordered[1:]:
  if g.is_constant():break
  g=g.gcd(a)
 if not g.is_constant():r={j:v//g for j,v in r.items()}
 lead=min(r.values(),key=len).leading_coefficient()
 if lead!=1:r={j:v/lead for j,v in r.items()}
 return r
def clear_entry(row,piv,j):
 if j not in row:return row
 a,b=row[j],piv[j]
 if a==b: ai=bi=a**0
 else:
  g=a.gcd(b); ai=a//g;bi=b//g
 out={k:bi*v for k,v in row.items() if k!=j}
 for k,v in piv.items():
  if k!=j:
   x=out.get(k,0)-ai*v
   if x:out[k]=x
   elif k in out:del out[k]
 return primitive(out)
def fraction_free_echelon(rows,limit=None):
 if not rows:return {}
 typ=next((type(x) for row in rows for x in row.values()),None)
 work=[cleared(r) for r in rows if r]
 piv={}
 while work:
  j=min(min(r) for r in work)
  idx=min((i for i in range(len(work)) if j in work[i]),key=lambda i:(len(work[i]),sum(len(v) for v in work[i].values())))
  p=work.pop(idx)
  work=[t for row in work if (t:=clear_entry(row,p,j))]
  piv[j]=p
  if limit==len(piv):break
 for j in sorted(piv,reverse=True):
  for k in sorted(piv[j]):
   if k>j and k in piv:piv[j]=clear_entry(piv[j],piv[k],k)
 return {j:{k:typ(v,row[j]) for k,v in row.items()} for j,row in piv.items()}

if __name__=='__main__':
    names=('b','c','v','w')
    field=rational_function_field(names)
    params=[field.from_sympy(v) for v in sp.symbols(names)]
    basis,image=field_paths(field_relations(*params),field.one,fraction_free_echelon)
    assert all([sum(len(basis[p,n,h]) for h in range(4)) for n in range(6)]==[1,4,10,20,35,56]
               for p in range(4))
    for h,expected in [(0,3),(1,2)]:
        equations,size=central_equations(basis,image,h)
        rank=len(fraction_free_echelon(equations))
        assert size-rank==expected
        print('Four independent parameters, character',h,'columns',size,'rank',rank,flush=True)
