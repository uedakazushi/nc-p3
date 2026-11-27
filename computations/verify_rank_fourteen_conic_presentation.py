"""All-degree regularity on the full period-two curve of horizontal 112 at opposite rank 14."""
from functools import cache
from itertools import combinations
import json
from pathlib import Path

from flint import nmod_mat
from gaussian_curve_field import ONE, convert, echelon, h, one, polynomial_keys, strip
from quartic_center_presentation import certify_presentation
from verify_klein_center_lift import field_paths,central_equations
from verify_parallel_curve_residues import residue_test
from verify_parallel_horizontal_moduli import window
from verify_parallel_invertible_rectangles import K,matrix,s
from verify_second_rectangle_low_fibers import b as bb,t as tt,x as xx,y as yy

FIELDS = Path(__file__).resolve().parents[1]/'research'/'rank_fourteen_conic_residues.json'


@cache
def conic_data():
    i = s.I
    A,B,C,D = h*h+2*h-2,h*h-2*h-2,h*h+2*i*h+2,h*h-2*i*h+2
    u,v,z = -A*C/(B*D),-A*D/(B*C),i*(h**4-4)/(4*h*h)
    def simp(q):
        return s.factor(q,extension=i)
    t = simp((u+v)/(1+u*v))
    b,x = simp(u*t),simp(u*(z+v)/(1+u*v))
    raw = window(tt,bb,(tt**2-bb)/(1-bb),xx,yy)
    W = {k:simp(f.xreplace({tt:t,bb:b,xx:x,yy:x})) for k,f in raw.items()}
    boundary = [u,v,u+v,1+u*v,z*z-1,z+v,z-u,u*z-1,v*z+1]
    def multiples(inputs):
        unit = ONE
        for f in inputs:
            if f:
                a = convert(f)
                for g in (a.norm(),a.d):
                    if g.degree()>0:
                        unit *= g//g.gcd(unit)
        return unit
    base = multiples(boundary)
    full = multiples(list(W.values())+boundary)
    assert all(strip(convert(f).d,base).degree()<=0 for f in W.values())
    return raw,W,base,full


def divisors(alternate):
    raw,W,base,full = conic_data()
    units = base if alternate else full
    shifts = (2,2,0,0) if alternate else (0,0,2,2)
    rel = []
    for j in K:
        rows = []
        for g in K:
            _,idx = matrix(raw,j,g).rref()
            M = matrix(W,j,g)[:,list(idx)]
            polys = [convert(s.factor(M.extract(J,range(M.cols)).det(),extension=s.I)).norm()
                     for J in combinations(K,M.cols)]
            gcd = next(q for q in polys if q)
            for q in polys:
                gcd = gcd.gcd(q)
            assert strip(gcd,units).degree()==0
            rows.extend((g^shifts[j]^shifts[(j+1)%4],
                         {(k^shifts[j],k^g^shifts[(j+1)%4]):convert(M[k,n]) for k in K if M[k,n]})
                        for n in range(M.cols))
        rel.append(rows)
    excluded = set()
    def span(rows):
        return echelon(rows,units,excluded)
    def pivot(rows):
        indices = []
        return echelon(rows,units,excluded,indices),indices
    basis,image = field_paths(rel,one,span)
    sections = [[] for _ in K]
    for g in (0,1):
        rows,n = central_equations(basis,image,g)
        P = span(rows)
        free = [i for i in range(n) if i not in P]
        assert len(free)==3-g
        for j in free:
            coeff = {j:one}|{k:-row[j] for k,row in P.items() if j in row}
            offset = 0
            for k in K:
                words = basis[k,4,g]
                sections[k].append((g,{w:coeff[offset+i] for i,w in enumerate(words) if offset+i in coeff}))
                offset += len(words)
    certify_presentation(rel,sections,0,[one*int(i==int(alternate)) for i in range(5)],
                         one,pivot,span,(1,4,10,20,30,36,34,20,0))
    assert len(excluded)==(171 if alternate else 174)
    return excluded|polynomial_keys(units,base)


def modular(h,i,p):
    A,B,C,D = h*h+2*h-2,h*h-2*h-2,h*h+2*i*h+2,h*h-2*i*h+2
    u,v = -A*C*pow(B*D,-1,p)%p,-A*D*pow(B*C,-1,p)%p
    z = i*(h**4-4)*pow(4*h*h,-1,p)%p
    assert u*v*(u+v)*(u*v+1)*(z*z-1)*(z+v)*(z-u)*(u*z-1)*(v*z+1)%p
    t,x = s.Rational(u+v,1+u*v),s.Rational(u*(z+v),1+u*v)
    def cv(f):
        n,d = s.fraction(f)
        return s.Integer(int(n)*pow(int(d),-1,p)%p)
    W = {k:cv(f) for k,f in window(t,u*t,v*t,x,x).items()}
    rules = []
    for j in K:
        rows = []
        for g in K:
            M = matrix(W,j,g)
            R,rank = nmod_mat([[int(x) for x in row] for row in M.tolist()],p).rref()
            assert rank==(2,2,1,1)[g]
            idx = [next(k for k,x in enumerate(row) if x) for row in R.tolist()[:rank]]
            rows.extend((g,{(k,k^g):int(M[k,n]) for k in K if M[k,n]}) for n in idx)
        rules.append(rows)
    return rules


def check():
    data = json.loads(FIELDS.read_text())
    assert divisors(False)&divisors(True)=={tuple(f) for f,_,_,_ in data}
    for f,p,r,i in data:
        poly = s.Poly.from_list(list(reversed(list(map(s.Rational,f)))),gens=h)
        assert poly.is_irreducible and poly.eval(r)%p==0 and poly.diff().eval(r)%p and (i*i+1)%p==0
        for j in (i,-i):
            residue_test(r,p,lambda value,prime:modular(value,j,prime))
    print('PASS: both full presentations and both embeddings at every one of the fourteen residue factors')


if __name__=='__main__':
    check()
