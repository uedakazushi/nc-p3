"""Sparse exact linear algebra over QQ(u,v) for the cubic-frame boundary.

No numerical sampling or polynomial interpolation is used for field ranks.
"""
from itertools import combinations_with_replacement,product
import sympy as s
from frobenius_tensor import digits

u,v=s.symbols('u v')
K=s.QQ.frac_field(u,v)
zero,one=K.zero,K.one


def tensor():
    monomials=list(combinations_with_replacement(range(4),3))
    lookup={m:i for i,m in enumerate(monomials)}
    rows=[]
    for a in range(4):
        for i in range(3):
            row=[0]*20
            row[lookup[tuple(sorted((i,i,a)))]]+=1
            row[lookup[tuple(sorted((i,3,a)))]]-=1
            rows.append(row)
        row=[0]*20
        for i,c in enumerate((1,u,v)):
            row[lookup[tuple(sorted((i,3,a)))]]+=c
        rows.append(row)
    equations=s.Matrix(rows)
    free=[(0,1,2),(1,3,3),(2,3,3),(3,3,3)]
    normalization=s.Matrix(4,20,lambda i,j:int(monomials[j]==free[i]))
    square=equations.col_join(normalization)
    assert s.factor(square.det())==-2*u*v
    basis=(square.inv()*s.zeros(16,4).col_join(s.eye(4))).applyfunc(s.cancel)
    assert (equations*basis).applyfunc(s.cancel)==s.zeros(16,4)
    assert normalization*basis==s.eye(4)
    result={}
    for encoded in range(256):
        a,b,j,c=digits(encoded)
        value=K.from_sympy(basis[lookup[tuple(sorted((a,b,c)))],j])
        if value:result[encoded]=value
    return result

def echelon(rows):
 pivots={}
 for row in rows:
  col={i:x for i,x in row.items() if x}
  while col:
   i=min(col)
   if i not in pivots:
    scale=col[i];pivots[i]={j:x/scale for j,x in col.items()};break
   scale=col[i]
   for j,x in pivots[i].items():
    y=col.get(j,zero)-scale*x
    if y:col[j]=y
    else:col.pop(j,None)
 return pivots

def nullspace(rows,size):
 piv=echelon(rows);answer=[]
 for free in range(size):
  if free in piv:continue
  z={free:one}
  for i,row in sorted(piv.items(),reverse=True):
   value=-sum((x*z.get(j,zero) for j,x in row.items() if j!=i),zero)
   if value:z[i]=value
  answer.append(z)
 return answer

def inserted(windows,start,n):
 for pos in range(n-1):
  for left in product(range(4),repeat=pos):
   for right in product(range(4),repeat=n-pos-2):
    for rel in windows[(start+pos)%4]:
     col={}
     for (a,b),value in rel.items():
      index=0
      for letter in left+(a,b)+right:index=4*index+letter
      col[index]=value
     yield col

def windows(w):
    result=[]
    for pos in range(4):
        columns={}
        for encoded,value in w.items():
            word=digits(encoded)
            a=4*word[pos]+word[(pos+1)%4]
            b=4*word[(pos+2)%4]+word[(pos+3)%4]
            columns.setdefault(b,{})[a]=value
        basis=echelon(columns.values())
        assert len(basis)==6
        result.append([{divmod(j,4):x for j,x in col.items()} for col in basis.values()])
    return result


def dual_windows(relations):
    result=[]
    for window in relations:
        rows=[{4*a+b:value for (a,b),value in row.items()} for row in window]
        kernel=nullspace(rows,16)
        assert len(kernel)==10
        result.append([{divmod(i,4):x for i,x in row.items()} for row in kernel])
    return result


def quotient(relations,start,degree):
    pivots=echelon(inserted(relations,start,degree))
    free=sorted(set(range(4**degree))-set(pivots))
    # Reuse the echelon rows to recover the orthogonal kernel.
    basis=nullspace(pivots.values(),4**degree)
    assert len(basis)==len(free)
    assert all(row.get(j,zero)==int(i==k)
               for i,row in enumerate(basis) for k,j in enumerate(free))
    return basis,free


def sparse_full_pivot_rank(columns,progress=None):
    """Exact column elimination, choosing low-complexity FLINT field pivots.

    Only the pivot order differs from echelon(); this routine returns a rank,
    not a nullspace or a monomial quotient frame.
    """
    remaining=[{i:x for i,x in column.items() if x} for column in columns]
    remaining=[column for column in remaining if column]
    rank=0
    while remaining:
        row_count={}
        for column in remaining:
            for row in column:row_count[row]=row_count.get(row,0)+1
        best=None;chosen=None
        for j,column in enumerate(remaining):
            for i,value in column.items():
                score=(len(value.n)+len(value.d),value.n.total_degree()+value.d.total_degree(),
                       (len(column)-1)*(row_count[i]-1),len(column))
                if best is None or score<best:
                    best,chosen=score,(j,i)
        j,i=chosen
        pivot_column=remaining.pop(j)
        pivot=pivot_column[i]
        normalized={row:x/pivot for row,x in pivot_column.items() if row!=i}
        updated=[]
        for column in remaining:
            scale=column.pop(i,zero)
            if scale:
                for row,value in normalized.items():
                    total=column.get(row,zero)-scale*value
                    if total:column[row]=total
                    else:column.pop(row,None)
            if column:updated.append(column)
        remaining=updated;rank+=1
        if progress is not None and rank%16==0:progress(rank)
    return rank
