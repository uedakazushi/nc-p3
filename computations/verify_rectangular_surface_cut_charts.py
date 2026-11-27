"""Regenerate the exceptional divisors of the four generic central cuts on surface 145."""
from functools import cache
from pathlib import Path
import json
import sys

from flint import fmpz_mpoly_ctx
from flint_rational_function import rational_function_field
from localized_unit_echelon import unit_echelon
from verify_klein_center_lift import field_paths,central_equations
from verify_parallel_rectangle_surface import K,U,a,matrix,s,window,y
from verify_rank_fourteen_central_cut import FOUR_QUARTIC_PROFILE,quotient_layers

CERTIFICATE = Path(__file__).resolve().parents[1]/'research'/'rectangular_surface_cut_factors.json'


def central_cut_data(w,field,units):
    rel = [[(h,{(j,j^h):field.from_sympy(s.factor(c[j])) for j in K if c[j]})
            for h in K for c in matrix(w,p,h).columnspace()] for p in K]
    return relation_cut_data(rel,field,units)


def relation_cut_data(rel,field,units,full_center=False,reducer=None):
    def reduce(rows):
        return unit_echelon(rows,units) if reducer is None else reducer(rows)
    basis,image = field_paths(rel,field.one,reduce)
    sections = [[] for _ in K]
    for h in (0,1):
        rows,width = central_equations(basis,image,h)
        piv = reduce(rows)
        free = [j for j in range(width) if j not in piv]
        assert len(free)==3-h
        for j in ((free[0],free[2]) if h==0 and not full_center else free):
            coefficients = {j:field.one}|{k:-row[j] for k,row in piv.items() if j in row}
            offset = 0
            for p in K:
                words = basis[p,4,h]
                sections[p].append((h,{w:coefficients[offset+i] for i,w in enumerate(words)
                                       if offset+i in coefficients}))
                offset += len(words)
            assert offset==width
    return field,units,rel,sections


def cut_data():
    field = rational_function_field(('a','y'))
    return central_cut_data(window(a,y,1),field,field.from_sympy(U*(a*y+1)).n)


def chart(p,data,profile=FOUR_QUARTIC_PROFILE):
    field,units,rel,sections = data
    factors = set()
    integer_ring = fmpz_mpoly_ctx.get(units.context().names())
    @cache
    def extract(expression):
        # Factors over Q are stored as primitive integer polynomials.
        from flint import fmpq_mpoly
        poly = fmpq_mpoly(expression,units.context())
        for f,_ in poly.factor()[1]:
            n = integer_ring.from_dict({k:int(c) for k,c in f.to_dict().items()})
            assert all(c.denominator==1 for c in f.to_dict().values())
            if n.leading_coefficient()<0:
                n = -n
            factors.add(str(n//int(n.content())))
    def record(f):
        extract(str(f))
    def project(n,h,rows,size):
        piv = unit_echelon(rows,units,record)
        free = [i for i in range(size) if i not in piv]
        ids = {j:i for i,j in enumerate(free)}
        return free,[{ids[j]:field.one} if j in ids else {ids[k]:-v for k,v in piv[j].items() if k!=j}
                     for j in range(size)]
    limit = len(profile)-1
    assert list(quotient_layers(rel,sections,p,limit,field.one,project,profile))==list(range(1,limit+1))
    return sorted(factors)


if __name__=='__main__':
    recorded = json.loads(CERTIFICATE.read_text())
    assert list(map(len,recorded))==[64,127,106,109]
    data = cut_data()
    for p in (list(map(int,sys.argv[1:])) or K):
        assert chart(p,data)==recorded[p]
        print('PASS: all central-cut divisors at start',p,flush=True)
