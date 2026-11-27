"""Eliminate the remaining zero-dimensional exceptions on the rectangular surface."""
from functools import cache
import json
import sys

from flint import fmpz_mpoly,fmpz_mpoly_ctx,fmpz_poly,nmod_poly
from sympy import primerange

from verify_rectangular_surface_cut_charts import CERTIFICATE
from verify_rectangular_surface_curve_strata import check_columns,good

RING = fmpz_mpoly_ctx.get(('a','y'))
PRIME = 1000003
CURVES = (
    {'a + 1','y + 1','a*y - 1','a*y - a + 1','a*y - y + 1','a*y^2 + 1'},
    {'a + 1','y + 1','a*y - 1','a*y - a + 1','a*y - y + 1','a*y^2 + 1',
     'a - y','a*y - 2*a + 1','a*y - 2*y + 1','a*y^2 - 1'},
)


def signature(f):
    if f.leading_coefficient()<0:
        f = -f
    return tuple((ij,int(n)) for ij,n in (f//int(f.content())).to_dict().items())


@cache
def polynomial(f):
    return RING.from_dict(dict(f))


@cache
def transform(f,swap,invert):
    poly = polynomial(f)
    da,dy = poly.degrees()
    out = {}
    for (i,j),n in f:
        if invert:
            i,j = da-i,dy-j
        out[(j,i) if swap else (i,j)] = n
    return signature(RING.from_dict(out))


def open_sets(start):
    charts = json.loads(CERTIFICATE.read_text())
    out = []
    for p,swap in ((start,False),(start+2,True)):
        fs = [signature(fmpz_mpoly(v,RING)) for v in charts[p]]
        out.extend({transform(f,swap,invert) for f in fs} for invert in (False,True))
    return sorted(out,key=len)


@cache
def projection(f,g):
    result = polynomial(f).resultant(polynomial(g),'y').to_dict()
    assert result and all(j==0 for i,j in result)
    degree = max(i for i,j in result)
    u = fmpz_poly([int(result.get((i,0),0)) for i in range(degree+1)])
    u = u//int(u.content())
    u = u//u.gcd(u.derivative()) if degree else fmpz_poly([1])
    if u.leading_coefficient()<0:
        u = -u
    for factor in (fmpz_poly([0,1]),fmpz_poly([-1,1])):
        if u%factor==0:
            u = u//factor
    return tuple(map(int,u.coeffs()))


@cache
def modular_polynomial(coeff):
    return nmod_poly(list(coeff),PRIME)


@cache
def common(left,right):
    if min(len(left),len(right))==1:
        return ()
    a,b = modular_polynomial(left),modular_polynomial(right)
    if a.degree()==len(left)-1 and b.degree()==len(right)-1 and a.gcd(b).degree()==0:
        return ()
    gcd = fmpz_poly(list(left)).gcd(fmpz_poly(list(right)))
    return tuple(map(int,gcd.coeffs())) if gcd.degree()>0 else ()


def abscissas(f,sets):
    pools = sorted(({projection(f,g) for g in group} for group in sets),key=len)
    keep = {a for a in pools[0] if len(a)>1}
    for pool in pools[1:]:
        keep = {h for a in keep for b in pool if (h:=common(a,b))}
    return {tuple(map(int,h.coeffs())) for a in keep for h,_ in fmpz_poly(list(a)).factor()[1]}


@cache
def fix_abscissa(f,x,p):
    coeff = [0]*(1+max(j for (i,j),n in f))
    for (i,j),n in f:
        coeff[j] = (coeff[j]+n*pow(x,i,p))%p
    return nmod_poly(coeff,p)


def boundary_data(f,A):
    # Work before reduction modulo p. The three deleted ordinate values are
    # distinct in characteristic zero here, since the abscissa -1 is handled
    # by its already certified full curve.
    equation = fmpz_poly(list(A))
    poly = polynomial(f)
    dx,dy = map(int,poly.degrees())
    terms = dict(f)
    for degree in range(dy,-1,-1):
        if fmpz_poly([terms.get((i,degree),0) for i in range(dx+1)])%equation:
            break
    else:
        raise AssertionError('A nonvertical irreducible curve contains a vertical line')
    multiplicities = []
    for test in ('zero','one','negative_reciprocal'):
        derivative = poly
        for order in range(degree+1):
            coeffs = derivative.to_dict()
            n = max((j for i,j in coeffs),default=0)
            terms1 = {}
            for (i,j),value in coeffs.items():
                if test=='zero' and j:
                    continue
                k = i+n-j if test=='negative_reciprocal' else i
                sign = -1 if test=='negative_reciprocal' and j%2 else 1
                terms1[k] = terms1.get(k,0)+sign*int(value)
            value = fmpz_poly([terms1.get(i,0) for i in range(max(terms1,default=0)+1)])
            if value%equation:
                multiplicities.append(order)
                break
            derivative = derivative.derivative('y')
        else:
            raise AssertionError('Multiplicity exceeds degree')
    return degree,multiplicities


@cache
def residue_abscissas(A,p):
    u = nmod_poly(list(A),p)
    if u.degree()!=len(A)-1:
        return []
    return [int(z) for z,e in u.roots() if e==1 and int(z) not in (0,1,p-1)]


@cache
def good_orbit(point,p):
    a,y = point
    for x,z in ((a,y),(y,a),(pow(a,-1,p),pow(y,-1,p)),(pow(y,-1,p),pow(a,-1,p))):
        if good((x,z),p):
            return True
    return False


def certify_fiber(f,A,swap,sets):
    if A==(1,1):
        return  # The minus curve and its half-turn have been certified.
    degree,multiplicities = boundary_data(f,A)
    for p in primerange(1009,5000):
        for x in residue_abscissas(A,p):
            possible = fix_abscissa(f,x,p)
            if possible.degree()!=degree:
                continue
            for z,n in zip((0,1,-pow(x,-1,p)%p),multiplicities):
                linear = nmod_poly([-z,1],p)
                if n:
                    if possible%(linear**n)!=0:
                        break
                    possible = possible//(linear**n)
                if not possible(z):
                    break
            else:
                for group in sets:
                    if possible.degree()==0:
                        break
                    product = nmod_poly([1],p)
                    for g in group:
                        product = (product*(fix_abscissa(g,x,p)%possible))%possible
                        if not product:
                            break
                    if product:
                        possible = possible.gcd(product)
                if possible.degree()==0:
                    return
                roots = possible.roots()
                if sum(e for _,e in roots)!=possible.degree():
                    continue
                if all(good_orbit((int(z),x) if swap else (x,int(z)),p) for z,e in roots):
                    return
    raise AssertionError(('Unverified number-field fiber',str(polynomial(f)),A))


def check(start):
    groups = open_sets(start)
    assert {str(polynomial(f)) for f in set.intersection(*groups)}==CURVES[start]
    count = total_degree = 0
    for f in sorted(groups[0],key=lambda z:(polynomial(z).total_degree(),str(polynomial(z)))):
        missing = [g for g in groups[1:] if f not in g]
        if not missing:
            continue
        swap = polynomial(f).degrees()[1]==0
        ff = transform(f,swap,False)
        sets = [{transform(g,swap,False) for g in group} for group in missing]
        for A in sorted(abscissas(ff,sets)):
            certify_fiber(ff,A,swap,sets)
            count += 1
            total_degree += len(A)-1
    assert (count,total_degree)==((721,2592),(1070,3938))[start]
    print('PASS: all residual fibers for start',start,'and its half-turn',flush=True)


if __name__=='__main__':
    check_columns()
    for start in (list(map(int,sys.argv[1:])) or (0,1)):
        check(start)
