"""Dimension bounds for endpoint projections, verified over all their integer ranks."""
from itertools import product


def rectangular_lower(d,k,l):
    x = 4-d
    first_pair = 4*x-(6-k)
    triple = 16*x-(6-k)*4-x*6+(4-l)
    return triple*l-first_pair*k


assert [rectangular_lower(1,3,l) for l in (2,3,4)] == [13,30,45]
assert min(rectangular_lower(2,k,l) for k,l in product((2,3,4),repeat=2)) == 4
assert 3*4*2-4*2-3*5 == 1
admissible = [(k,h,j) for k,h,j in product(range(2,7),range(2,7),range(17))
              if 16-2*(6-k)-2*(6-h)+j <= 1]
assert admissible == [(2,2,0),(2,2,1)]
print("PASS: rectangular bounds, one-dimensional endpoint cut and two-line overlap obstruction")
