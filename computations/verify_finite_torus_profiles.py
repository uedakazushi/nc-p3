"""Enumerate additive coincidences and local complete-intersection lengths at fixed points."""
from collections import Counter
from functools import reduce
from itertools import combinations, product
from math import gcd, prod

import sympy as sp



def primitive(v):
    divisor = reduce(gcd, v)
    if not divisor:
        return tuple(v)
    first = next(x for x in v if x)
    return tuple((1 if first > 0 else -1)*x//divisor for x in v)


def determinant3(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) + a[1]*(b[2]*c[0]-b[0]*c[2])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


def candidate_gaps():
    row = [(0,0,0,0), (1,0,0,0), (1,0,1,0), (1,1,1,0)]
    col = [(0,0,0,0), (1,0,0,0), (1,0,0,1), (1,1,0,1)]
    cells = [tuple(x+y for x,y in zip(a,b)) for a,b in product(row,col)]
    planes = {primitive(tuple(x-y for x,y in zip(a,b))) for a,b in combinations(cells,2)}
    planes = {v for v in planes if min(v) < 0 < max(v)}
    lines = set()
    for triple in combinations(planes,3):
        v = primitive(tuple((-1)**i*determinant3(*(r[:i]+r[i+1:] for r in triple)) for i in range(4)))
        if min(v)>0:
            lines.add(v)
    assert len(planes) == 18 and len(lines) == 31
    return lines


def fixed_records(a, b):
    counts = Counter(x+y for x,y in product(a,b))
    active = [z for z,n in counts.items() if n>1]
    for multiplicities in product(*(range(counts[z]) for z in active)):
        if sum(multiplicities) != 6:
            continue
        rel = dict(zip(active,multiplicities))
        fixed = [(i,j) for i,x in enumerate(a) for j,y in enumerate(b) if not rel.get(x+y)]
        if min(len({i for i,j in fixed}),len({j for i,j in fixed})) != len(fixed):
            continue
        lengths = []
        for i,j in fixed:
            weight = a[i]+b[j]
            numerator = prod((h-weight)**v for h,v in rel.items())
            denominator = prod(x-a[i] for x in a if x!=a[i])*prod(y-b[j] for y in b if y!=b[j])
            length = sp.Rational(numerator,denominator)
            if not length.is_Integer or length <= 0:
                break
            lengths.append(int(length))
        if len(lengths) == len(fixed) and sum(lengths) == 20:
            yield tuple(sorted(h for h,v in rel.items() for _ in range(v))),tuple(zip(fixed,lengths))


def check_edges():
    result = {}
    for d,e,s,t in candidate_gaps():
        a, b = (0,d,d+s,d+s+e),(0,d,d+t,d+t+e)
        for rec in fixed_records(a,b):
            assert (d,e,s,t) not in result
            result[d,e,s,t] = rec
    fixed = (((0,0),10),((3,3),10))
    assert result == {(1,1,1,1): ((1,2,3,3,4,5),fixed),
                      (1,1,1,2): ((1,2,3,4,5,6),fixed),
                      (1,1,2,1): ((1,2,3,4,5,6),fixed)}


def check_alternation():
    t = sp.Symbol('t')
    a = ((0,1,2,3),(0,1,3,4))*2
    char = [sum(t**x for x in v) for v in a]
    rel = sum(t**i for i in range(1,7))
    # The coefficient of t^2 in the first power sum is negative.
    first = (prod(char) - sum(rel*char[(i+2)%4]*char[(i+3)%4] for i in range(4))
             + t**7*sum(x*x.subs(t,1/t) for x in char) + 2*rel**2 - 4*t**7)
    assert first.expand().coeff(t,2) == -2
    local = (1-t)**2*(1-t**2)**2*(1-t**3)**2/((1-t)**4*(1-t**2)**2)
    assert sp.cancel(local - (1+t+t**2)**2) == 0


if __name__ == '__main__':
    check_edges()
    check_alternation()
    print('Finite fixed schemes: four simple weights have two length-ten sites; a double middle weight has 9+2+9.')
