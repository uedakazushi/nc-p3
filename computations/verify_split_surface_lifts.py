"""Positive-weight Schur equations over the generic Bell surface Q(u,v)."""
from itertools import product
from subprocess import run

import sympy as sp

from sympy.polys.matrices import DomainMatrix

from verify_isolated_split_point_scheme import local_pencil
from verify_split_inverse_windows import BITS, PAIRS, bell, borders, reciprocal, transport

WORDS = tuple(product(range(4),repeat=4))
BINARY_PAIRS = tuple(product(range(4),repeat=2))
U,V = sp.symbols('u v')


def split_tensor(h):
    p = bell([sp.Integer(1),-U*V,U,V])
    z = sp.diag(1,-1)
    maps = ([sp.eye(2)]*4,[h*z,h,h*z,h])
    targets = [transport(reciprocal(p,i),maps[i],PAIRS[i]) for i in range(2)]
    pivot = next(q for q in BITS if targets[0][q] != 0)
    scale = sp.cancel(targets[0][pivot]/targets[1][pivot])
    assert all(sp.cancel(targets[0][q]-scale*targets[1][q]) == 0 for q in BITS)
    colors = (0,1,0,1)
    def word(q,pattern):
        return tuple(b+2*c for b,c in zip(q,pattern))
    w = {word(q,colors):x for q,x in p.items()}
    w.update({word(q,[1-c for c in colors]):x for q,x in targets[0].items()})
    for i,q in product(range(2),BITS):
        values = borders(maps[i],PAIRS[i],i,q)
        for k in (0,1):
            pattern = [c ^ (int(j not in {i,i+1}) ^ k) for j,c in enumerate(colors)]
            w[word(q,pattern)] = values[k]*(scale if i == 1 and k == 0 else 1)
    return w


PROCEDURES = """
LIB "linalg.lib";
ring r=(0,u,v),dummy,dp;
proc linearConstraints(matrix l,matrix r,intvec rows,intvec cols)
{
  int i,j,k;
  matrix a[100][size(rows)];
  for(k=1;k<=size(rows);k++){
    for(i=1;i<=10;i++){
      for(j=1;j<=10;j++){a[(i-1)*10+j,k]=l[i,rows[k]]*r[cols[k],j];}
    }
  }
  return(a);
}
proc shaped(matrix a,int k,intvec rows,intvec cols)
{
  matrix b[16][16]; int i;
  for(i=1;i<=size(rows);i++){b[rows[i],cols[i]]=a[i,k];}
  return(b);
}
proc stack(matrix a,matrix b){return(transpose(concat(transpose(a),transpose(b))));}
proc vec(matrix a)
{
  matrix b[100][1]; int i,j;
  for(i=1;i<=10;i++){for(j=1;j<=10;j++){b[(i-1)*10+j,1]=a[i,j];}}
  return(b);
}
"""

CHECK_SLICE = """
matrix ker=syz(d3);
matrix span=g;
matrix z[64][4];
int i,j,n;
intvec idx=1..64;
matrix best;
for(i=1;i<=ncols(ker);i++){
  best=concat(span,submat(ker,idx,intvec(i)));
  if(rank(best)>rank(span)){n++;z[1..64,n]=ker[1..64,i];span=best;}
}
print("DIMS "+string(ncols(ker))+" "+string(rank(d4))+" "+string(rank(g))+" "+string(n));
matrix rhs[200][10]; n=0;
matrix x0,x1,y0,y1,s0,s1,rt;
for(i=1;i<=4;i++){for(j=i;j<=4;j++){
  n++;
  x0=shaped(z,i,p03,q03); x1=shaped(z,i,p13,q13);
  y0=shaped(z,j,p03,q03); y1=shaped(z,j,p13,q13);
  s0=x0*j0*y0; s1=x1*j1*y1;
  if(i!=j){s0=s0+y0*j0*x0;s1=s1+y1*j1*x1;}
  rt=stack(vec(l0*s0*r0),vec(l1*s1*r1));
  rhs[1..200,n]=rt[1..200,1];
}}
print("AUGMENTED "+string(rank(concat(d4,rhs))));
quit;
"""


def singular_input(w):
    code = [PROCEDURES]
    def matrix(name,rows):
        data = ','.join(str(sp.factor(a)).replace('**','^') for row in rows for a in row)
        code.append(f'matrix {name}[{len(rows)}][{len(rows[0])}]={data};')
    def vec(name,indices):
        code.append('intvec '+name+'='+','.join(str(i+1) for i in indices)+';')
    for i in range(2):
        def rotate(q):
            return q[-i:] + q[:-i] if i else q
        rows = [[w.get(rotate(a+b),0) for b in BINARY_PAIRS] for a in BINARY_PAIRS]
        # The rational coordinate merely chooses a nonzero pivot minor.
        sample = sp.Matrix(rows).subs({U:2,V:3})
        row_pivots,col_pivots = sample.T.rref()[1],sample.rref()[1]
        assert len(row_pivots) == len(col_pivots) == 6
        matrix(f'm{i}',rows)
        vec(f'a{i}',row_pivots)
        vec(f'b{i}',col_pivots)
        code.append(f'''matrix l{i}=transpose(syz(transpose(m{i})));matrix r{i}=syz(m{i});
matrix e{i}=inverse(submat(m{i},a{i},b{i}));matrix j{i}[16][16];
int ii{i},kk{i};
for(ii{i}=1;ii{i}<=6;ii{i}++){{for(kk{i}=1;kk{i}<=6;kk{i}++){{
 j{i}[b{i}[ii{i}],a{i}[kk{i}]]=e{i}[ii{i},kk{i}];
}}}}''')
        for n in (3,4):
            words = [q[i:]+q[:i] for q in WORDS if sum(j<2 for j in q) == n]
            vec(f'p{i}{n}',[4*q[0]+q[1] for q in words])
            vec(f'q{i}{n}',[4*q[2]+q[3] for q in words])
    for n in (3,4):
        code.append(f'matrix d{n}=stack(linearConstraints(l0,r0,p0{n},q0{n}),'
                    f'linearConstraints(l1,r1,p1{n},q1{n}));')
    words = [q for q in WORDS if sum(j<2 for j in q) == 3]
    matrix('g',[[w.get(q[:i]+(t,)+q[i+1:],0) if q[i] == n else 0
                 for i,n,t in product(range(4),range(2),range(2,4))] for q in words])
    return '\n'.join(code+[CHECK_SLICE])


def verify_graph(w):
    coeffs = {q:sp.QQ.from_sympy(sp.sympify(v).subs({U:2,V:3})) for q,v in w.items()}
    def block(i,j,k,l):
        return DomainMatrix.from_list([[coeffs.get((a,b,i,j),sp.QQ.zero)
                                        for b in range(l,l+2)] for a in range(k,k+2)],sp.QQ)
    def relation(k,l):
        return next(m for i,j in product(range(l,l+2),repeat=2)
                    if not (m:=block(i,j,k,k)).is_zero_matrix)
    nn,tt = relation(0,2),relation(2,0)
    nt,tn = ([block(i,j,k,l) for i,j in product(range(2),range(2,4))]
             for k,l in ((0,2),(2,0)))
    assert local_pencil(sp.QQ,nn,tt,nt,tn) == local_pencil(sp.QQ,tt,nn,tn,nt) == 0


def verify():
    x,z = sp.Matrix([[0,1],[1,0]]),sp.diag(1,-1)
    for h in (sp.eye(2),x,z,x*z):
        w = split_tensor(h)
        verify_graph(w)
        result = run(['Singular','-q'],input=singular_input(w),text=True,
                     capture_output=True,check=True)
        assert '?' not in result.stdout and 'DIMS 20 13 16 4' in result.stdout, result.stdout
        assert 'AUGMENTED 13' in result.stdout, result.stdout
    print("Generic positive slices of all four completions have dimensions 4 and 3 and zero obstruction.")


if __name__ == '__main__':
    verify()
