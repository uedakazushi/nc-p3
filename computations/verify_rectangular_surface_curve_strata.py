"""The residual rational curves of the central-cut certificates on the rectangle surface."""
from functools import cache
import sys

from flint import fmpz_poly
from sympy import primerange
from verify_parallel_curve_residues import residue_test
from verify_parallel_double_zero_curve import chart
from verify_parallel_rectangle_surface import K,a,matrix,s,window,y


def columns(a,y):
    b = a*y
    d,e = a+y-b,a+y-1
    return (
        ((a,d,1,1),(b,y,1,0),(a,1,1,0),(e,y,1,1),(1,1,1,1),(1,1,1,0)),
        ((a,e,1,1),(1,y,1,0),(y,d,1,1),(b,a,0,1),(1,1,1,1),(0,1,1,1)),
        ((a,b,1,0),(d,y,1,1),(a,e,1,1),(1,y,1,0),(1,1,1,0),(1,1,1,1)),
        ((a,1,1,0),(e,y,1,1),(y,b,0,1),(d,a,1,1),(1,0,1,1),(1,1,1,1)))


def check_columns():
    w = window(a,y,1)
    for p,rows in enumerate(columns(a,y)):
        expected = [list(c) for h in K for c in matrix(w,p,h).columnspace()]
        assert list(map(list,rows))==expected
        for i in (0,2):
            assert abs(s.Matrix([rows[i][2:],rows[i+1][2:]]).det())==1


def modular_rules(point,prime):
    return [[(h,{(j,j^h):int(value%prime) for j,value in enumerate(col) if value%prime})
             for h,col in zip((0,0,1,1,2,3),rs)] for rs in columns(*point)]


@cache
def good(point,prime):
    try:
        residue_test(point,prime,modular_rules)
        return True
    except AssertionError:
        return False


# Entries give a as a function of y, the parameter involution, and the number
# of factors before/after reflection and the sum of the remaining degrees.
CURVES = {
    'minus':(lambda v:-s.S.One,'reciprocal',(1911,67,136)),
    'diagonal':(lambda v:v,'reciprocal',(2234,48,126)),
    'unit-product':(lambda v:1/v,'reciprocal',(1988,47,88)),
    'graph':(lambda v:1/(1-v),1,(2341,68,200)),
    'graph2':(lambda v:1/(2-v),2,(2425,58,124)),
    'square':(lambda v:1/v**2,'reciprocal',(2427,46,124)),
    'negative-square':(lambda v:-1/v**2,'reciprocal',(2429,39,98)),
}


def run_curve(label):
    first,move,counts = CURVES[label]
    def system(v,field):
        # Take the column space after specializing: on a boundary curve the
        # pivot choice can change while the relation subspace has constant rank.
        w = window(first(v),v,1)
        return [[(h,{(i,i^h):field.from_sympy(s.factor(c[i])) for i in K if c[i]})
                 for h in K for c in matrix(w,p,h).columnspace()] for p in K]
    factors = chart(system)
    def reflected(f):
        poly = (fmpz_poly(list(reversed(f))) if move=='reciprocal'
                else fmpz_poly(list(f))(fmpz_poly([move,-1])))
        return tuple(map(int,(poly if poly.leading_coefficient()>0 else -poly).coeffs()))
    residual = [f for f in factors if f[0] and reflected(f) in factors]
    assert (len(factors),len(residual),sum(len(f)-1 for f in residual))==counts
    v = s.Symbol('v')
    u = first(v)
    num,den = s.fraction(s.factor(u*v*(u-1)*(v-1)*(u*v+1)))
    omitted = {tuple(map(int,reversed(s.Poly(f,v).all_coeffs())))
               for z in (num,den) for f,_ in s.factor_list(z)[1]}
    # Clear any additional parameter pole, even if the displayed product cancels.
    omitted.update(tuple(map(int,reversed(s.Poly(f,v).all_coeffs())))
                   for f,_ in s.factor_list(s.denom(s.cancel(u)))[1])
    for f in sorted(set(residual)-omitted):
        poly = fmpz_poly(list(f))
        passed = False
        for p in primerange(1009,2200):
            if not f[-1]%p:
                continue
            for z in range(2,p):
                if int(poly(z))%p:
                    continue
                n,d = s.fraction(first(s.Integer(z)))
                if int(d)%p:
                    x = int(n)*pow(int(d),-1,p)%p
                    if x not in (0,1) and (x*z+1)%p and good((x,z),p):
                        passed = True
                        break
            if passed:
                break
        assert passed,(label,f)
    print('PASS: full residual curve',label,flush=True)


if __name__=='__main__':
    check_columns()
    for label in (sys.argv[1:] or CURVES):
        run_curve(label)
