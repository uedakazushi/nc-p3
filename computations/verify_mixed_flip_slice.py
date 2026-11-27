"""All torus tangent jumps and the full scheme germ at the mixed component intersection."""
import sympy as sp
from sympy.polys.matrices import DomainMatrix

from verify_split_simple_surface import V, EXTRA, D0
from verify_split_inverse_windows import bell
from verify_isolated_split_component import PRIME, WORDS, residue, nmod_mat
from verify_split_spectral_line import I, multiply, reconstruct
from verify_split_spectral_line import A,line_tensor,torus
from itertools import product
from functools import reduce
from collections import defaultdict
from subprocess import run

def label(i,j):return (j^int(i>=2)) if j<2 else (j^(i%2))
def degree(q):return tuple(sum(x<2 and label(i,x)==a for i,x in enumerate(q))-1 for a in (0,1))
def kernel_blocks():
 w=line_tensor(sp.I,multiply(torus(V,sp.Integer(1)),A))
 field=sp.QQ_I.frac_field(V);tensor={q:field.from_sympy(sp.sympify(v)) for q,v in w.items()};ks=[]
 for i in range(4):
   pairs=sorted(product(range(4),repeat=2),key=lambda q:label(i,q[0])<=label((i+1)%4,q[1]))
   def val(q):return tensor.get(q[-i:]+q[:-i] if i else q,field.zero)
   m,p=DomainMatrix.from_list([[val(x+y) for x in pairs] for y in pairs],field).rref();assert p==tuple(range(6));r=m.to_list()
   ks.append([{q:(-r[j][k] if j<6 else field.one) for j,q in enumerate(pairs) if (r[j][k] if j<6 else j==k)} for k in range(6,16)])
 rows=defaultdict(list);cols=defaultdict(list)
 for q in product(range(4),repeat=4):cols[degree(q)].append(q)
 for i in (0,1):
  for a,b in product(ks[i],ks[i+2]):
   row={}
   for x,u in a.items():
    for y,v in b.items():
     q=x+y;row[q[-i:]+q[:-i] if i else q]=u*v
   h=degree(next(iter(row)));assert all(degree(q)==h for q in row)
   den=reduce(lambda a,b:a.lcm(b),(v.denom for v in row.values()))
   rows[h].append({q:(v.numer*(den//v.denom)).as_expr() for q,v in row.items()})
 return [(h,cols[h],rows[h]) for h in sorted(cols)]
def smith_characters():
 script=['LIB "jacobson.lib";ring r=(0,c),v,dp;minpoly=c^2+1;']
 for h,cols,rows in kernel_blocks():
    size=len(cols)
    vals=','.join(str(r.get(q,0)).replace('**','^').replace('I','c') for r in rows for q in cols)
    script.append(f'print("CHAR {h[0]} {h[1]}");matrix m[{len(rows)}][{size}]={vals};')
    if size==1:
        script.append('poly d=m[1,1];int j;for(j=2;j<=nrows(m);j++){d=gcd(d,m[j,1]);};'
                      'if(deg(d)>0 || d==0){print("ERROR");};print("ZERO 0");kill d,j,m;')
    else:
        script.append('matrix h=smith(transpose(m));int j;int z=0;for(j=1;j<=ncols(m);j++)'
                      '{if(j>ncols(h)){z++;}else{if(h[j,j]==0){z++;}else{poly f=h[j,j];'
                      'while(f%v==0){f=f/v;};if(deg(f)>0){print("FACTOR "+string(f/leadcoef(f)));};'
                      'kill f;};};};print("ZERO "+string(z));kill j,z,h,m;')
 out=run(['Singular','-q'],input='\n'.join(script+['quit;']),text=True,capture_output=True,check=True).stdout
 assert '?' not in out and 'ERROR' not in out
 ranks, factors = {},defaultdict(list)
 for line in out.splitlines():
    if line.startswith('CHAR'):
        _,i,j=line.split()
        key=int(i),int(j)
    elif line.startswith('ZERO'):
        ranks[key]=int(line.split()[1])
    elif line.startswith('FACTOR'):
        factors[key].append(line.split()[1])
 assert sum(ranks.values())==60 and ranks[0,0]==20
 assert {k for k,v in factors.items() if v==['v2-1']}=={(-1,1),(1,-1)}
 roots={(-1,0),(0,-1),(-1,2),(2,-1),(0,2),(1,1),(2,0)}
 assert {k for k,v in factors.items() if sorted(v)==sorted(['v+(c)','v+(-c)']) or v==['v2+1']}==roots
 assert len(factors)==9


def planes():
 x,y,a=sp.symbols('x y a')
 p=bell((sp.Integer(1),sp.Integer(1),sp.I,sp.I))
 p[0,0,1,1]=x;p[1,0,0,1]=sp.I*x
 p[1,1,0,0]=y;p[0,1,1,0]=sp.I*y
 extra=multiply(EXTRA,D0)
 plane=reconstruct(p,extra)
 line=reconstruct(bell((sp.Integer(1),sp.Integer(1),sp.I,sp.I)),multiply(torus(sp.Integer(1),a),extra))
 def sample(w,der=None,vals=None):
    vv={q:sp.sympify(v) for q,v in w.items()}
    return {q:residue((z.diff(der) if der is not None else z).subs(vals or {})) for q,z in vv.items()}
 base=sample(line,vals={a:1})
 plane_base=sample(plane,vals={x:0,y:0})
 assert all((base.get(q,0)-plane_base.get(q,0))%PRIME==0 for q in WORDS)
 return base,sample(line,a,{a:1}),sample(plane,x,{x:0,y:0}),sample(plane,y,{x:0,y:0})
def obstruction():
 tensors=planes();eqs=[]
 for i in (0,1):
    m,z,x,y=[nmod_mat(16,16,[w.get(q[-i:]+q[:-i] if i else q,0) for q in WORDS],PRIME) for w in tensors]
    left,n=m.transpose().nullspace();right,b=m.nullspace();assert n==b==10
    def pivots(c):
        r,n=c.rref();return [next(k for k in range(c.ncols()) if int(r[j,k])) for j in range(n)]
    cols=pivots(m);rows=pivots(m.transpose());inv=nmod_mat([[int(m[j,k]) for k in cols] for j in rows],PRIME).inv()
    g=nmod_mat(16,16,PRIME)
    for j,k in product(range(6),repeat=2):g[cols[j],rows[k]]=inv[j,k]
    assert m*g*m==m
    rhs=[left.transpose()*(z*g*d+d*g*z)*right for d in (x,y)]
    eqs.extend([[int(left[4*q[i]+q[(i+1)%4],j]*right[4*q[(i+2)%4]+q[(i+3)%4],k]) for q in WORDS]+[int(v[j,k]) for v in rhs] for j,k in product(range(10),repeat=2)])
 assert [nmod_mat([row[:k] for row in eqs],PRIME).rank() for k in (256,257,258)]==[194,195,196]
 w,*directions=tensors
 orbit=[[w.get(q[:i]+(b,)+q[i+1:],0) if q[i]==a else 0 for i,a,b in product(range(4),repeat=3)]
        +[v.get(q,0) for v in directions] for q in WORDS]
 assert [nmod_mat([row[:k] for row in orbit],PRIME).rank() for k in range(64,68)]==[59,60,61,62]


if __name__=='__main__':
    smith_characters()
    obstruction()
    print('All slice tangent weights and both independent mixed quadratic obstructions verified.')
