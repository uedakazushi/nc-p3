"""Rank certificates over the entire invertible ordinary Bell surface."""
from itertools import product
from subprocess import run

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from verify_split_surface_lifts import U, V, WORDS, split_tensor

FIELD = sp.QQ.frac_field(U, V)
PAIR = tuple(product(range(4), repeat=2))
UNIT = U*V*(U**2-1)*(V**2-1)
PHASE = (sp.I, sp.I, -sp.I, -sp.I)


def polynomial_row(row):
    row = [sp.factor(a) for a in row]
    denominator = sp.factor(sp.lcm([sp.denom(a) for a in row]))
    for factor, _ in sp.factor_list(denominator)[1]:
        assert sp.denom(sp.cancel(UNIT/factor)).is_number
    return tuple(sp.factor(a*denominator) for a in row)


def equations(tensor):
    kernels = []
    for i in range(4):
        def val(q):
            return tensor.get(q[-i:]+q[:-i] if i else q, FIELD.zero)
        pivots = [(a,b+2) for a,b in product(range(2), repeat=2)]
        for offset in (0, 2):
            pivots.append(next(p for p in PAIR if all(offset <= a < offset+2 for a in p)
                               and any(val(p+t) for t in PAIR)))
        pairs = pivots + [q for q in PAIR if q not in pivots]
        m = DomainMatrix.from_list([[val(x+y) for x in pairs] for y in PAIR], FIELD)
        reduced, cols = m.rref()
        assert cols == tuple(range(6))
        values = reduced.to_list()
        kernel = [{q: (-values[j][k] if j < 6 else FIELD.one) for j,q in enumerate(pairs)
                   if (values[j][k] if j < 6 else j == k)} for k in range(6,16)]
        for row in kernel:
            polynomial_row([FIELD.to_sympy(a) for a in row.values()])
        kernels.append(kernel)
    rows = {2: [], 3: []}
    for i in (0, 1):
        for a,b in product(kernels[i], kernels[i+2]):
            row = {}
            for x,l in a.items():
                for y,r in b.items():
                    q = x+y
                    row[q[-i:]+q[:-i] if i else q] = l*r
            n = sum(j < 2 for j in next(iter(row)))
            assert all(sum(j < 2 for j in q) == n for q in row)
            if n in rows:
                rows[n].append(row)
    assert [len(rows[n]) for n in (2, 3)] == [68, 48]
    return rows


# Eliminate only pivots invertible on the indicated open set. Finish by testing
# the ideal of maximal required minors in the ring localized at loc.
ELIMINATE = r"""
LIB "elim.lib";
ring r=0,(u,v,z),dp;
proc rankBound(matrix a,int needed,poly loc)
{
 int i,j,k,l,x,y; poly p,g,d;
 while(needed>0){
    x=0;
    for(i=1;i<=nrows(a) && x==0;i++){for(j=1;j<=ncols(a) && x==0;j++){
       if(a[i,j]!=0 && deg(a[i,j])==0){x=i;y=j;}
    }};
    for(i=1;i<=nrows(a) && x==0;i++){for(j=1;j<=ncols(a) && x==0;j++){
       p=a[i,j];if(p!=0){
          d=gcd(p,loc);while(deg(d)>0){p=p div d;d=gcd(p,loc);};
          if(deg(p)==0){x=i;y=j;}
       }
    }};
    if(x==0){break;};
    matrix b[nrows(a)-1][ncols(a)-1];k=0;
    for(i=1;i<=nrows(a);i++){if(i!=x){k++;l=0;g=0;
       for(j=1;j<=ncols(a);j++){if(j!=y){l++;
          b[k,l]=a[x,y]*a[i,j]-a[x,j]*a[i,y];g=gcd(g,b[k,l]);
       }};
       if(g!=0){d=gcd(g,loc);while(deg(d)>0){
          for(l=1;l<=ncols(b);l++){b[k,l]=b[k,l] div d;};
          g=g div d;d=gcd(g,loc);
       };};
    }};
    a=b;kill b;needed--;
 }
 if(needed>0){
    ideal ii=minor(a,needed);ii=ii,1-z*loc;ii=std(ii);
    if(reduce(1,ii)!=0){print("FAILED");return();}
 }
 print("VERIFIED");
}
poly base=u*v*(u^2-1)*(v^2-1);
poly smooth=base*(u+v)*(u*v-1)*(u^2+1);
"""


def verify():
    w = {q: FIELD.from_sympy(sp.sympify(a)) for q,a in split_tensor(sp.eye(2)).items() if a}
    rows = equations(w)
    script = [ELIMINATE]

    def certify(matrix, rank, locus):
        vals = ','.join(str(a).replace('**','^') for row in matrix for a in row)
        script.append(f'matrix a[{len(matrix)}][{len(matrix[0])}]={vals};'
                      f'rankBound(a,{rank},{locus});kill a;')

    for n, parity, flip in product((2,3), (0,1), (1,-1)):
        words = [q for q in WORDS if sum(j < 2 for j in q) == n and sum(q)%2 == parity and q[0]%2 == 0]
        def phase(q):
            return sp.prod(PHASE[i] for i,j in enumerate(q) if j < 2)/(flip*(sp.I if n==3 else 1))
        block = []
        for row in rows[n]:
            vals = polynomial_row([
                FIELD.to_sympy(row.get(q,FIELD.zero))
                + phase(q)*FIELD.to_sympy(row.get(tuple(j^1 for j in q),FIELD.zero))
                for q in words])
            if any(vals) and vals not in block and tuple(-a for a in vals) not in block:
                block.append(vals)
        if n == 2:
            rank = 18 if (parity, flip) == (0, 1) else 16
            # Two opens cover both the simple spectrum and the nonline no-flip
            # collision u=v (where u^4 != 1).
            certify(block, rank, 'base*(u^2-v^2)*(u^2*v^2-1)')
            certify(block, rank, 'smooth')
        else:
            certify(block, 11, 'base')
            shear = []
            for q in words:
                other = tuple(j^1 for j in q)
                def coef(p, i):
                    key = p[:i]+(2+parity,)+p[i+1:]
                    return FIELD.to_sympy(w.get(key, FIELD.zero)) if p[i] == 0 else 0
                shear.append(polynomial_row([sp.factor(coef(q,i)+phase(q)*coef(other,i))
                                             for i in range(4)]))
            certify(shear, 4, 'base')
    result = run(['Singular','-q'],input='\n'.join(script+['quit;']),text=True,
                 capture_output=True,check=True).stdout
    assert '?' not in result and 'FAILED' not in result and result.count('VERIFIED') == 16, result


if __name__ == '__main__':
    verify()
    print('All tangent and shear rank bounds verified on the simple and no-flip collision loci.')
