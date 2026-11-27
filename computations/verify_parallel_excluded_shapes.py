"""Full coefficient exclusion for the asymmetric parallel thin rectangles."""
from functools import lru_cache
from itertools import product

import sympy as s
from sympy.polys.matrices import DomainMatrix

from flint_rational_function import rational_function_field
from verify_klein_generic_quartics import fraction_free_echelon
from verify_klein_mixed_positive import matrix
from verify_klein_parallel_flags import relaxed_fourth
from verify_klein_sparse_transverse import torus_forms
from verify_klein_three_lines import K
from verify_klein_two_doublet_supports import candidates


@lru_cache(None)
def factors(poly):
    return frozenset(f for f,_ in s.factor_list(poly)[1] if not f.is_Symbol)


def cells(equations,units=frozenset(),sub=None):
    if sub is None:
        sub = {}
    pool = {factors(f)-units for f in equations if f!=0}
    if frozenset() in pool:
        return
    pool = {fs for fs in pool if not any(g<fs for g in pool)}
    if not pool:
        yield sub
        return
    multiple = [fs for fs in pool if len(fs)>1]
    if multiple:
        fs = min(multiple,key=lambda g:(len(g),tuple(sorted(map(str,g)))))
        rest = [s.prod(g) for g in pool if g!=fs]
        for f in sorted(fs,key=str):
            yield from cells(rest+[f],units,sub)
            units = units|{f}
        return
    equations = [next(iter(fs)) for fs in pool]
    options = []
    for f in equations:
        for x in f.free_symbols:
            poly = s.Poly(f,x)
            if poly.degree()!=1:
                continue
            A,B = s.factor(poly.nth(1)),s.factor(poly.nth(0))
            val = s.factor(-B/A)
            key = (bool((factors(A)-units) and (factors(B)-units)),s.count_ops(val),
                   sum(x in g.free_symbols for g in equations),str(x),str(val))
            options.append((key,f,x,A,B,val))
    assert options
    key,f,x,A,B,val = min(options,key=lambda row:row[0])
    rest = [g for g in equations if g!=f]
    if key[0]:
        yield from cells(rest+[A,B],units,sub)
    replacement = {x:val}
    def numerator(poly):
        return s.fraction(s.factor(poly.xreplace(replacement)))[0]
    new_units = [numerator(g) for g in units|{A,B}]
    if 0 in new_units:
        return
    new_sub = {j:s.factor(a.xreplace(replacement)) for j,a in sub.items()}
    new_sub[x] = val
    yield from cells([numerator(g) for g in rest],frozenset().union(*(factors(g) for g in new_units)),new_sub)


def zero_fourth(w):
    index = {(i,j,k,i^j^k):n for n,(i,j,k) in enumerate(product(K,repeat=3))}
    rows = []
    for pos,h in product(range(3),K):
        M = matrix(w,pos,h)
        pivots = DomainMatrix.from_Matrix(M).to_field().rref()[1]
        if len(pivots)<(2 if h<2 else 1):
            # A generic upper bound for the whole block, hence also for every fiber.
            return 'rank'
        assert len(pivots)==(2 if h<2 else 1)
        for k in pivots:
            col = M[:,k]
            for j in K:
                outer = (j,j^h)
                r = [0]*64
                for i in K:
                    r[index[outer[:pos]+(i,i^h)+outer[pos:]]]=col[i]
                rows.append(r)
    return 64-DomainMatrix.from_Matrix(s.Matrix(rows)).to_field().rank()


def exceptional_chart(eq,y):
    a,b,c,d,e = y
    C = d*(a*b+1)/a
    E2 = (a*b+1)*(b*d+1)/a**2
    rest = [s.factor(f.subs(c,C).subs(e**2,E2)) for f in eq]
    F = (a-1)*(a*b+1)*d+a**2
    assert all(g==0 or not factors(s.fraction(s.factor(g/F))[0])-{a*b+1} for g in rest)
    assert any(s.factor(g*a/(b*d)-F)==0 for g in rest)
    D = -a**2/((a-1)*(a*b+1))
    B = (a-1)*(1-a**2*e**2)/a
    assert s.factor(E2.subs(d,D).subs(b,B)-e**2)==0
    result = {b:B,d:s.factor(D.subs(b,B)),c:s.factor(C.subs(d,D).subs(b,B))}
    assert all(s.factor(f.xreplace(result))==0 for f in eq)
    return result


DIRECT = {73:1,86:1,87:1,91:0,92:0,93:0,97:0,98:0,100:0,102:1,103:0,104:0}
RUN = (*range(67,76),86,87,*range(91,94),*range(97,112))


def path_dimension(w,start,length,final_color):
    sy = sorted(set().union(*(v.free_symbols for v in w.values())),key=str)
    field = rational_function_field(tuple(map(str,sy)) or ('v',))
    relations = {}
    for p,b in product(K,K):
        M = matrix(w,p,b)
        echelon,pivots = DomainMatrix.from_Matrix(M.T).to_field().rref()
        if len(pivots)<1+(b<2):
            return None
        assert len(pivots)==1+(b<2)
        E = echelon.to_Matrix()
        relations[p,b] = [[field.from_sympy(s.factor(v)) for v in E.row(j)] for j in range(len(pivots))]
    dimensions = [[1,0,0,0]]
    previous = {}
    for n in range(1,length+1):
        current,arrows = [],{}
        for h in (K if n<length else (final_color,)):
            offsets = []
            size = 0
            for j in K:
                offsets.append(size)
                size += dimensions[-1][h^j]
            rows = []
            if n>1:
                for b in K:
                    g = h^b
                    for col in relations[(start+n-2)%4,b]:
                        for a in range(dimensions[-2][g]):
                            row = {}
                            for i,coefficient in enumerate(col):
                                if coefficient:
                                    j = i^b
                                    row.update({offsets[j]+k:coefficient*v for k,v in previous[g,i][a].items()})
                            if row:
                                rows.append(row)
            piv = fraction_free_echelon(rows)
            free = [i for i in range(size) if i not in piv]
            current.append(len(free))
            position = {j:k for k,j in enumerate(free)}
            projection = [{position[i]:field.one} if i in position else
                          {position[k]:-v for k,v in piv[i].items() if k!=i} for i in range(size)]
            for j in K:
                arrows[h^j,j] = projection[offsets[j]:offsets[j]+dimensions[-1][h^j]]
        previous = arrows
        dimensions.append(current)
    return dimensions[-1][0]


LATE = {79:(0,7,0,30),83:(3,7,1,30),123:(0,7,0,30),128:(0,7,0,30),
        130:(1,5,1,14),131:(2,6,3,22),137:(3,7,1,30),141:(3,7,1,30)}


IMPOSSIBLE = (80,81,84,119,120,121,122,125,126,127,132,133,135,136,138,140)


def exceptional_planes():
    p,q,a,b,V,T = s.symbols('p q a b V T')
    M = s.Matrix([[a,0,0,V],[0,-p*a-b,-V,T],[0,-V,b,0],[V,T,0,-a-q*b]])
    det = s.factor(M.det())
    assert s.factor(det.subs(V,0)-a*b*((p*a+b)*(a+q*b)-T**2))==0
    assert s.factor(det.subs({p:0,q:0,T:0})-(a**2+V**2)*(b**2+V**2))==0
    assert s.factor(det.subs({q:1/p,V:0})-a*b*((p*a+b)**2/p-T**2))==0


def main():
    exceptional_planes()
    for i in IMPOSSIBLE:
        assert all(eq==[1] for _,_,eq in torus_forms(candidates(1)[i],1))
    for i in (*RUN,*LATE):
        for w,y,eq in torus_forms(candidates(1)[i],1):
            if eq==[1]:
                continue
            if i in DIRECT:
                assert relaxed_fourth(w,DIRECT[i])>9
                continue
            charts = [exceptional_chart(eq,y)] if i==107 else cells(eq)
            for sub in charts:
                assert all(s.factor(f.xreplace(sub))==0 for f in eq)
                W = {k:s.factor(a.xreplace(sub)) for k,a in w.items()}
                if i in LATE:
                    p,n,g,bound = LATE[i]
                    dim = path_dimension(W,p,n,g)
                    assert dim is None or dim>bound, (i,dim)
                else:
                    dim = zero_fourth(W)
                    assert dim=='rank' or dim>9, (i,dim)
    print('PASS: every torus phase and all residual coefficient strata on 29 quartic and eight later exclusions, 16 further empty tori, and both six-node planes')


if __name__=='__main__':
    main()
