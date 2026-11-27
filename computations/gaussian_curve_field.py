"""Exact arithmetic and local pivot certificates over Q(i)(h), using real polynomial norms."""
from collections import Counter
from functools import cache

from flint import fmpq, fmpq_poly
import sympy as s

h = s.Symbol('h')
ONE = fmpq_poly([1])
ZERO = fmpq_poly([])


class Gaussian:
    __slots__ = ('a','b','d')

    def __init__(self,a,b=ZERO,d=ONE):
        if not d:
            raise ZeroDivisionError
        if not a and not b:
            self.a,self.b,self.d = ZERO,ZERO,ONE
            return
        if not d.is_one():
            common = a.gcd(b).gcd(d)
            if common.degree()>0:
                a,b,d = a//common,b//common,d//common
            lead = d.leading_coefficient()
            if lead!=1:
                a,b,d = a/lead,b/lead,d/lead
        self.a,self.b,self.d = a,b,d

    def __bool__(self):
        return bool(self.a or self.b)

    def __add__(self,other):
        other = coerce(other)
        if not other:
            return self
        if not self:
            return other
        if self.d==other.d:
            return Gaussian(self.a+other.a,self.b+other.b,self.d)
        g = self.d.gcd(other.d)
        left,right = self.d//g,other.d//g
        return Gaussian(self.a*right+other.a*left,self.b*right+other.b*left,self.d*right)

    __radd__ = __add__

    def __neg__(self):
        return Gaussian(-self.a,-self.b,self.d)

    def __sub__(self,other):
        return self+-coerce(other)

    def __mul__(self,other):
        other = coerce(other)
        return Gaussian(self.a*other.a-self.b*other.b,self.a*other.b+self.b*other.a,self.d*other.d)

    __rmul__ = __mul__

    def __truediv__(self,other):
        other = coerce(other)
        return Gaussian((self.a*other.a+self.b*other.b)*other.d,
                        (self.b*other.a-self.a*other.b)*other.d,self.d*other.norm())

    def __eq__(self,other):
        other = coerce(other)
        return (self.a,self.b,self.d)==(other.a,other.b,other.d)

    def norm(self):
        """Numerator norm; its zeros contain the numerator zeros for either embedding of i."""
        return self.a*self.a+self.b*self.b

    def weight(self):
        return max(self.a.degree(),self.b.degree(),self.d.degree())


def coerce(value):
    if isinstance(value,Gaussian):
        return value
    return Gaussian(fmpq_poly([fmpq(value)]))


one,zero = coerce(1),coerce(0)


@cache
def convert(expr):
    def polynomial(value):
        coeff = [s.expand_complex(c).as_real_imag() for c in reversed(s.Poly(value,h,extension=s.I).all_coeffs())]
        def part(j):
            return fmpq_poly([fmpq(int(t[j].p),int(t[j].q)) for t in coeff])
        return Gaussian(part(0),part(1))
    num,den = s.fraction(s.cancel(expr,extension=s.I))
    return polynomial(num)/polynomial(den)


def strip(f,units):
    while f.degree()>0 and (common:=f.gcd(units)).degree()>0:
        f = f//common
    return f


def polynomial_keys(f,units):
    f = strip(f,units)
    if f.degree()<=0:
        return set()
    return {tuple(map(str,p.coeffs())) for p,_ in f.factor()[1]}


def echelon(rows,units,excluded,indices=None):
    def clear(row,pivot,j):
        if j in row:
            c = -row.pop(j)
            for k,a in pivot.items():
                if k!=j:
                    val = row.get(k,zero)+c*a
                    if val:
                        row[k] = val
                    else:
                        row.pop(k,None)
        return row

    work = [(i,{j:x for j,x in row.items() if x}) for i,row in enumerate(rows) if any(row.values())]
    pivots = {}
    while work:
        counts = Counter(j for _,row in work for j in row)
        _,i,j = min(((int(x.weight()>0),(len(row)-1)*(counts[j]-1),x.weight(),len(row)),i,j)
                    for i,(_,row) in enumerate(work) for j,x in row.items())
        index,row = work.pop(i)
        if indices is not None:
            indices.append(index)
        x = row[j]
        excluded.update(polynomial_keys(x.norm(),units),polynomial_keys(x.d,units))
        pivot = {k:y/x for k,y in row.items()}
        pivots[j] = pivot
        work = [(k,row) for k,row in work if clear(row,pivot,j)]
    for j in reversed(pivots):
        for k in list(pivots[j]):
            if k!=j and k in pivots:
                clear(pivots[j],pivots[k],k)
    return pivots
