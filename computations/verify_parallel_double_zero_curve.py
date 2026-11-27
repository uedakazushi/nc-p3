"""The entire unmarked double-zero curve, its index shift, and all degree certificates."""
from functools import cache
from itertools import combinations, permutations

from flint import fmpq_poly, fmpz_poly

from flint_rational_function import rational_function_field
from verify_klein_center_lift import central_equations, field_paths
from verify_klein_generic_quartics import fraction_free_echelon
from verify_klein_three_lines import character
from verify_parallel_third_rectangle import K, WORDS, matrix, product, s, smith, third_rectangle
from verify_rank_fourteen_central_cut import quotient_layers
from verify_transverse_rank_twelve_isomorphism import alternating_gradient


def identifications():
    v = s.Symbol('v',nonzero=True)
    w = third_rectangle(v,2/v,2)
    keys = sorted(k for k,f in w.items() if f)
    gauge = s.Matrix([[int(k[p]==i) for p,i in product(K,repeat=2)] for k in keys])
    diag,left,_ = smith(gauge)
    rows = left.tolist()[sum(bool(a) for a in diag.diagonal()):]
    for param,rotation,perms in [
        (-v,0,[tuple(i^h for i in K) for h in (3,2,2,3)]),
        (2/v,2,[(1,0,2,3)]*4),
    ]:
        target = {k[rotation:]+k[:rotation]:z.subs(v,param) for k,z in w.items() if z}
        factors = [s.factor(w[k]/target[tuple(pi[i] for i,pi in zip(k,perms))]) for k in keys]
        assert all(character(row,factors)==1 for row in rows)
    w = {k:w.get(k,0) for k in WORDS}
    g = alternating_gradient(w)
    norm = s.factor(sum(w[k]*g[k] for k in WORDS)/4)
    assert norm==-16
    L = [[(matrix(w,p,h)*matrix(g,p,h).T/norm).applyfunc(s.factor) for h in K] for p in K]
    vals = [[s.factor(sum((m*m)[i,i] for m in row)) for i in K] for row in L]
    assert s.factor(sum(vals[0])+s.Rational(157,16))==0
    assert s.factor(sum(vals[1])-(v**2+4/v**2+316)/16)==0
    partial = sum(a*a for a in vals[0])
    assert s.factor(partial-partial.subs(v,2/v)-s.Rational(55,1024)*(v**2-4/v**2))==0
    for p in K:
        right = [s.factor(sum((m*m)[i^h,i^h] for h,m in enumerate(L[(p-1)%4]))) for i in K]
        for i,j in combinations(K,2):
            diffs = [s.fraction(s.factor(a[i]-a[j]))[0] for a in (right,vals[p])]
            assert s.gcd(*diffs)==1
    support = set(keys)
    rotated = {k[1:]+k[:1] for k in keys}
    for p,q,r in product(permutations(K),repeat=3):
        last = {}
        for i,j,k,l in keys:
            z = p[i]^q[j]^r[k]
            if l in last and last[l]!=z:
                break
            last[l]=z
        else:
            assert len(set(last.values()))<4 or {(p[i],q[j],r[k],last[l]) for i,j,k,l in support}!=rotated
    expected = [-(2*v**2-5*v+4)/v,(2*v**2+5*v+4)/v,-1,-1]
    for h in K:
        m = s.Matrix([[w[i,j,i^h,j^h] for j in K] for i in K])
        assert s.factor(m.det()-expected[h])==0


def path_system(v,field):
    w = third_rectangle(v,2/v,2)
    return [[(h,{(j,j^h):field.from_sympy(s.factor(c[j])) for j in K if c[j]})
             for h in K for c in matrix(w,p,h).columnspace()] for p in K]


def chart():
    field = rational_function_field(('v',))
    windows = path_system(s.Symbol('v'),field)
    excluded = set()
    @cache
    def record_primitive(coefficients):
        for f,_ in fmpz_poly(list(coefficients)).factor()[1]:
            excluded.add(tuple(map(int,f.coeffs())))
    def record(f):
        if f.is_constant():
            return
        a = f.to_dict()
        z = fmpq_poly([a.get((j,),0) for j in range(max(k[0] for k in a)+1)]).numer()
        record_primitive(tuple(map(int,(z//z.content()).coeffs())))
    def reduce(rows):
        return fraction_free_echelon(rows,record=record)
    for win in windows:
        for _,r in win:
            for v in r.values():
                record(v.d)
    basis,image = field_paths(windows,field.one,reduce)
    sections = [[] for _ in K]
    for h,weights in ((0,((1,3,5),(2,7,11))),(1,((1,0),(0,1)))):
        rows,width = central_equations(basis,image,h)
        piv = reduce(rows)
        free = [i for i in range(width) if i not in piv]
        assert len(free)==3-h
        for choice in weights:
            values = {j:c*field.one for j,c in zip(free,choice) if c}
            for k,row in piv.items():
                z = -sum((values[j]*row[j] for j in free if j in values and j in row),field.zero)
                if z:
                    values[k]=z
            offset = 0
            for p in K:
                words = basis[p,4,h]
                sections[p].append((h,{w:values[offset+j] for j,w in enumerate(words) if offset+j in values}))
                offset += len(words)
            assert offset==width
    def project(n,h,rows,size):
        piv = reduce(rows)
        free = [j for j in range(size) if j not in piv]
        ids = {j:i for i,j in enumerate(free)}
        out = [{ids[j]:field.one} if j in ids else {ids[k]:-v for k,v in piv[j].items() if j!=k} for j in range(size)]
        return free,out
    for p in K:
        assert list(quotient_layers(windows,sections,p,13,field.one,project))==list(range(1,14))
    return excluded


def orbit_filter(excluded):
    def key(coeffs):
        f = fmpz_poly(coeffs)
        sign = 1 if f.leading_coefficient()>0 else -1
        return tuple(map(int,(f//(sign*f.content())).coeffs()))
    out = []
    for f in sorted(excluded):
        reciprocal = [v*2**(len(f)-1-i) for i,v in enumerate(reversed(f))]
        if f[0] and all(key([v*sign**i for i,v in enumerate(a)]) in excluded
                        for sign in (-1,1) for a in (f,reciprocal)):
            out.append(f)
    return out


if __name__=='__main__':
    from verify_parallel_curve_residues import certify
    identifications()
    all_factors = chart()
    exceptional = orbit_filter(all_factors)
    assert (len(all_factors),len(exceptional),sum(len(f)-1 for f in exceptional))==(1885,40,64)
    certify(exceptional)
    print('PASS: all double-zero fibers, their unmarked coordinate, half-shift and every central specialization')
