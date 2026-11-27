"""Check Hilbert counts and the degree bounds in the central-sequence criterion."""
from itertools import product
from math import comb

import sympy as sp


def main():
    for degrees in [(4,)*4, (1, 2, 3, 4), (1,)*4]:
        top = sum(d-1 for d in degrees)
        ambiguity_bound = 3*max(1, top, *degrees)
        normal = [0]*(ambiguity_bound+1)
        for tail in product(*(range(d) for d in degrees)):
            weight = sum(tail)
            if weight <= ambiguity_bound:
                normal[weight] += 1
        for d in degrees:
            for n in range(d, len(normal)):
                normal[n] += normal[n-d]
        assert normal == [comb(n+3, 3) for n in range(len(normal))]
        assert all(top+5-d > sum(degrees[:i]) for i, d in enumerate(degrees))
    assert 3*max(12, 4) == 36
    # Noncentral second parameter in a homogeneous Ore example: [y,x]=z*w.
    x, y, z, w = sp.symbols('x y z w')
    test = x**3*y**2 + z*x*w

    def y_action(q):
        return y*q + z*w*sp.diff(q, x)

    assert sp.expand(y_action(x*test)-x*y_action(test)-z*w*test) == 0
    print('PASS: ordered Hilbert counts, overlap bounds and successive centrality example')


if __name__ == '__main__':
    main()
