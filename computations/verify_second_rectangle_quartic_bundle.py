"""A single unit-pivot basis of quartic centers over all large horizontal 112 supports."""
from itertools import combinations,product
from flint_rational_function import rational_function_field
from localized_unit_echelon import unit_echelon
from verify_klein_center_lift import field_paths,central_equations
from verify_parallel_horizontal_moduli import window
from verify_parallel_invertible_rectangles import K,matrix,s


def chart(t,b,c,x,y):
    raw = window(t,b,c,x,y)
    variables = tuple(sorted(set.union(*(e.free_symbols for e in raw.values())),key=str))
    field = rational_function_field(tuple(map(str,variables)))
    E = (b+c)*x*y-b*c*(x+y)+b*b*(c-1)
    required = (t,b,c,b-1,b+c,E,x,y,x-b,y-b,x-1,y-1,c*x+b*(1-c),c*y+b*(1-c))
    units = s.prod({f for e in required if e!=0 for g in s.fraction(s.factor(e)) for f,_ in s.factor_list(g)[1]})
    U = field.from_sympy(units).n
    shifts = (0,0,2,2)
    rel = []
    for p in K:
        q = (p+1)%4
        rows = []
        for h in K:
            cols = matrix(raw,p,h).columnspace()
            M = s.Matrix.hstack(*cols)
            minors = [s.factor(M.extract(I,range(M.cols)).det()) for I in combinations(K,M.cols)]
            assert any(d and all(s.div(units,f,*variables)[1]==0 for a in s.fraction(d)
                                 for f,_ in s.factor_list(a)[1]) for d in minors)
            rows.extend((h^shifts[p]^shifts[q],
                         {(i^shifts[p],i^h^shifts[q]):field.from_sympy(s.factor(col[i])) for i in K if col[i]})
                        for col in cols)
        rel.append(rows)
    def reduce(rows):
        return unit_echelon(rows,U)
    basis,image = field_paths(rel,field.one,reduce)
    assert all(sum(len(basis[p,n,h]) for h in K)==(1,4,10,20,35,56)[n] for p,n in product(K,range(6)))
    for h,d in enumerate((3,2,0,0)):
        rows,size = central_equations(basis,image,h)
        assert size-len(reduce(rows))==d


if __name__=='__main__':
    b,t,x,y = s.symbols('b t x y')
    for sign in (s.S.One,-s.S.One):
        chart(sign,s.S.One,s.S.One,x,y)
    chart(t,b,(t*t-b)/(1-b),x,y)
    print('PASS: all path and character-resolved quartic-center fibers, using only common boundary units')
