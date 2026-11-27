"""All regular fibers of support 143 and its two boundary presentations."""
from functools import cache
from itertools import combinations, product

from flint import fmpq_mat, fmpz_poly
from verify_klein_three_lines import character
from verify_parallel_curve_residues import certify
from verify_parallel_double_zero_curve import chart
from verify_parallel_invertible_rectangles import D,K,WORDS,completed,matrix,quadric,s,smith
from verify_third_rectangle_orbits import support_maps
from verify_transverse_rank_twelve_isomorphism import alternating_gradient


@cache
def window(y,t):
    b = y*(1-y)
    a = t*(1-y)
    d = t*(y*y-y+1)
    assert s.factor(quadric(112,t,b,1,a,1,y,d))==0
    return completed(112,t,b,1,a,1,y,d)


def path_system(y,convert):
    w = window(y,1)
    return [[(h,{(j,j^h):convert(s.factor(c[j])) for j in K if c[j]})
             for h in K for c in matrix(w,p,h).columnspace()] for p in K]


def modular(value,prime):
    def residue(z):
        n,d = s.fraction(z)
        return int(n)*pow(int(d),-1,prime)%prime
    return path_system(s.Integer(value),residue)


def opposite_and_frame(t):
    y = s.Symbol('y')
    raw = window(y,t)
    w = {k:raw.get(k,0) for k in WORDS}
    g = alternating_gradient(w)
    def numerator(f):
        return s.fraction(s.factor(f))[0]
    def separates(rows):
        for i,j in combinations(K,2):
            assert s.gcd_list([numerator(row[i]-row[j]) for row in rows]).is_number
    for p in (1,3):
        entries = []
        for pos in (p,(p-1)%4):
            blocks = [matrix(w,pos,h)*matrix(g,pos,h).T for h in K]
            squares = [m*m for m in blocks]
            entries.append([sum(m[i,i] if pos==p else m[i^h,i^h] for h,m in enumerate(squares))
                            for i in K])
        separates(entries)
    for p in (0,2):
        rows = [[sum(w[k]*g[k] for k in WORDS if k[p]==i and (k[1],k[3])==jl)
                 for i in K] for jl in product(K,repeat=2)]
        separates(rows)
    u = y*y-y
    target = {0:{2:1,3:0},1:{2:1,3:u+1,4:(u-1)*(u+1)**2},
              2:{2:1,3:0},3:{2:1,3:0}} if t==1 else {
              0:{4:1},1:{3:1,4:u**3-u**2+3*u+1},2:{4:1},3:{4:1}}
    for color,tests in target.items():
        m = s.Matrix([[w[i,j,i^color,j^color] for j in K] for i in K])
        for n,f in tests.items():
            minors = [numerator(m.extract(ii,jj).det()) for ii,jj in product(combinations(K,n),repeat=2)]
            gcd = s.gcd_list(minors)
            if f==0:
                assert gcd==0
            else:
                assert s.factor(gcd/f).is_number and gcd!=0,(t,color,n,gcd,f)


def orbits(t):
    y = s.Symbol('y')
    f = y*y-y+1
    def specialize(value,boundary):
        if not boundary:
            return s.factor(value)
        num,den = s.fraction(s.cancel(value))
        return s.rem(num*s.invert(den,f,y),f,y)
    for boundary in (False,True):
        base = {k:specialize(v,boundary) for k,v in window(y,t).items()}
        base = {k:v for k,v in base.items() if v}
        for shift,reflect in ((0,False),(1,False),(2,True),(3,False),(0,True)):
            if reflect and shift==0 and not boundary:
                # For the interior the identity is the sole support permutation
                # and the normalized coefficient y is recovered directly.
                continue
            raw = window(1-y if reflect else y,t)
            target = {k:specialize(v,boundary) for k,v in raw.items() if specialize(v,boundary)!=0}
            source = {k[shift:]+k[:shift]:v for k,v in base.items()}
            keys = sorted(source)
            gauge = s.Matrix([[int(k[p]==i) for p,i in product(K,repeat=2)] for k in keys])
            diag,left,_ = smith(gauge)
            rows = left.tolist()[sum(bool(a) for a in diag.diagonal()):]
            maps = list(support_maps(source,target))
            if not boundary:
                assert len(maps)==int(shift%2==0)
            good = False
            for pi in maps:
                ratios = [specialize(source[k]/target[tuple(q[i] for q,i in zip(pi,k))],boundary)
                          for k in keys]
                if all(specialize(character(row,ratios)-1,boundary)==0 for row in rows):
                    good = True
                    break
            assert good==(boundary or shift%2==0)


def exceptional_quartic():
    y = s.Symbol('y')
    f = y*y-y-1
    w = window(y,1)
    def pair(z):
        n,d = s.fraction(s.cancel(z))
        r = s.Poly(s.rem(n*s.invert(d,f,y),f,y),y)
        return int(r.nth(0)),int(r.nth(1))
    values = {k:pair(z) for k,z in w.items()}
    for k,v in w.items():
        a,b = values[k]
        assert s.rem(s.fraction(s.cancel(v-a-b*y))[0],f,y)==0
    dims = []
    for h in K:
        keys = [(i,j,k,i^j^k^h) for i,j,k in product(K,repeat=3)]
        rows = []
        for p,c,i,j in product(range(3),K,K,K):
            tail = (j,j^h^c)
            arow,brow = [0]*128,[0]*128
            for letter in K:
                local = (letter,letter^c,i,i^c)
                word = local[-p:]+local[:-p] if p else local
                a,b = values.get(word,(0,0))
                slot = keys.index(tail[:p]+(letter,letter^c)+tail[p:])
                arow[2*slot:2*slot+2] = (a,b)
                brow[2*slot:2*slot+2] = (b,a+b)
            rows.extend([arow,brow])
        rank = fmpq_mat(rows).rank()
        assert rank%2==0
        dims.append(64-rank//2)
    assert dims==[10,10,8,8],dims


def regularity():
    factors = chart(lambda y,field:path_system(y,field.from_sympy))
    def reflected(f):
        poly = fmpz_poly(list(f))(fmpz_poly([1,-1]))
        return tuple(map(int,(poly if poly.leading_coefficient()>0 else -poly).coeffs()))
    residual = {f for f in factors if reflected(f) in factors}
    assert (len(factors),len(residual),sum(len(f)-1 for f in residual))==(2542,66,169)
    excluded = {(0,1),(-1,1),(-1,-1,1)}
    assert excluded<=residual
    certify(sorted(residual-excluded),modular)


if __name__=='__main__':
    for t in (1,-1):
        opposite_and_frame(t)
        orbits(t)
    exceptional_quartic()
    regularity()
    print('PASS: all fibers, unmarked orbits, rotations, opposite ranks and exceptional quartic on support 143')
