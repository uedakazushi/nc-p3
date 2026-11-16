"""Exact QQ(u,v) arithmetic using python-flint polynomial gcd.

Elements store a reduced numerator and monic denominator. No floating point,
finite-field evaluations, or interpolation are used. This is an optional
accelerator for the boundary verifier, not a different mathematical field.
"""
from flint import fmpq_mpoly_ctx,fmpq
import sympy as s
from functools import lru_cache
ctx=fmpq_mpoly_ctx.get(('u','v'))
unit=ctx.constant(1)

class Element:
 __slots__=('n','d')
 def __init__(self,n,d=unit):
  if not d:raise ZeroDivisionError
  if not n:self.n=ctx.constant(0);self.d=unit;return
  if d!=unit:
   g=n.gcd(d);n=n//g;d=d//g
   lead=d.leading_coefficient()
   if lead!=1:n=n/lead;d=d/lead
  self.n,self.d=n,d
 @staticmethod
 def coerce(other):
  if isinstance(other,Element):return other
  return Element(ctx.constant(fmpq(other)))
 def __bool__(self):return bool(self.n)
 def __eq__(self,other):
  other=self.coerce(other);return self.n==other.n and self.d==other.d
 def __add__(self,other):
  other=self.coerce(other)
  if not other:return self
  if not self:return other
  return Element(self.n*other.d+other.n*self.d,self.d*other.d)
 __radd__=__add__
 def __neg__(self):return Element(-self.n,self.d)
 def __sub__(self,other):return self+-self.coerce(other)
 def __rsub__(self,other):return self.coerce(other)+-self
 def __mul__(self,other):
  other=self.coerce(other)
  if not self or not other:return zero
  return Element(self.n*other.n,self.d*other.d)
 __rmul__=__mul__
 def __truediv__(self,other):
  other=self.coerce(other);return Element(self.n*other.d,self.d*other.n)
 def __rtruediv__(self,other):return self.coerce(other)/self
 def __repr__(self):return f'({self.n})/({self.d})'
 def as_expr(self):
  u,v=s.symbols('u v')
  def poly(p):return sum((s.Rational(str(c))*u**a*v**b for (a,b),c in p.to_dict().items()),s.S.Zero)
  return poly(self.n)/poly(self.d)
zero,one=Element(ctx.constant(0)),Element(unit)
class Field:
 zero,one=zero,one
 @staticmethod
 @lru_cache(None)
 def from_sympy(expr):
  u,v=s.symbols('u v');n,d=s.fraction(s.cancel(expr))
  def poly(x):return ctx.from_dict({mon:fmpq(int(c.p),int(c.q)) for mon,c in s.Poly(x,u,v).terms()})
  return Element(poly(n),poly(d))
K=Field()
