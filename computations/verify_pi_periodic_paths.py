"""Exact checks of path separation and of two periodic blocks joined through four phases."""
from math import prod

import sympy as sp


def separate():
    # Eight consecutive blocks of the distinct affine points (1,t,t^2,t^3).
    values = []
    for target in range(8):
        others = [j for j in range(8) if j != target]
        values.append([prod((8*row+i)-(8*j+i) for i,j in enumerate(others))
                       for row in range(8)])
    matrix = sp.Matrix(values)
    assert matrix.is_diagonal() and matrix.det() != 0


def projected_bridge():
    p,q,length = 4,8,5
    first = list(range(p))
    mid = list(range(p,p+length-1))
    last = list(range(p+length-1,p+length-1+q))
    size = p+q+length-1
    phase = [i%4 for i in first]+[(i+1)%4 for i in range(len(mid))]
    phase += [(length+i)%4 for i in range(q)]
    letters = [[sp.zeros(size) for _ in range(2)] for _ in range(4)]
    def arrow(i,j,scale,value):
        assert phase[j] == (phase[i]+1)%4
        letters[phase[i]][0][i,j] = scale
        letters[phase[i]][1][i,j] = scale*value
    for block,offset in ((first,10),(last,50)):
        for j,i in enumerate(block):
            arrow(i,block[(j+1)%len(block)],2 if i == last[-1] else 1,offset+j)
    chain = [first[0]]+mid+[last[0]]
    for j,(start,end) in enumerate(zip(chain,chain[1:])):
        arrow(start,end,1,20+j)
    x = sum((v[0] for v in letters), sp.zeros(size))
    a = x**8
    unit = sp.eye(size)
    assert a*(a-unit)*(a-2*unit) == sp.zeros(size)
    left,right = -a*(a-2*unit),a*(a-unit)/2
    transient = unit-left-right
    projectors = [transient*sp.diag(*[int(j==i) for j in phase]) for i in range(4)]
    assert all(e*e==e for e in [left,right]+projectors)
    assert all(projectors[i]*projectors[j]==sp.zeros(size) for i in range(4) for j in range(i))
    walk = left
    for i,start in enumerate(chain[:-1]):
        ordinary,other = letters[phase[start]]
        if i == 0:
            walk *= other-10*ordinary
        elif i == length-1:
            walk *= other-(50+q-1)*ordinary
        else:
            walk *= ordinary
        if i+1 < length:
            walk *= projectors[phase[chain[i+1]]]
    walk *= right
    assert walk[first[0],last[0]] == (20-10)*(20+length-1-(50+q-1))


if __name__ == '__main__':
    separate()
    projected_bridge()
    print("Eight point-block evaluations separate; the projected 4-to-8-cycle bridge is nonzero.")
