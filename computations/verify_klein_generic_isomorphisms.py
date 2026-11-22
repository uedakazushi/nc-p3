"""Check the projective configurations and the rational coordinates of the second sheet."""
from itertools import combinations, permutations
from collections import Counter
import sympy as sp
from verify_klein_two_doublets import cross_blocks, gram, modular_points, other_completion

P = 12619


def normalize(v,prime=P):
    k = pow(next(x for x in v if x%prime),-1,prime)
    return tuple(x*k%prime for x in v)


def kernel(rows,prime=P):
    for inds in combinations(range(rows.rows),3):
        m = rows[list(inds),:]
        v = [(-1)**i*int(m[:,[j for j in range(4) if j!=i]].det())%prime for i in range(4)]
        if any(v):
            assert all(int(x)%prime == 0 for x in rows*sp.Matrix(v))
            return normalize(v,prime)
    raise ValueError('Expected a one-dimensional point kernel')


def framed(points,prime=P):
    planes = set()
    for a,b,c in combinations(points,3):
        coeff = []
        for i in range(4):
            j,k,l = [n for n in range(4) if n!=i]
            det = (a[j]*(b[k]*c[l]-b[l]*c[k])-a[k]*(b[j]*c[l]-b[l]*c[j])
                   +a[l]*(b[j]*c[k]-b[k]*c[j]))
            coeff.append((-1)**i*det%prime)
        assert any(coeff)
        inds = tuple(i for i,q in enumerate(points) if sum(x*y for x,y in zip(coeff,q))%prime == 0)
        if len(inds)>3:
            planes.add(inds)
    assert len(planes)==48 and all(len(pl)==4 for pl in planes)
    incidences = [sum(i in pl for pl in planes) for i in range(20)]
    assert Counter(incidences) == {8:4, 10:16}
    basis = sp.Matrix.hstack(*(sp.Matrix(q) for q,n in zip(points,incidences) if n==8))
    inv = basis.inv_mod(prime)
    return [normalize(tuple(map(int,inv*sp.Matrix(q))),prime) for q in points]


def number_of_maps(source,target,prime=P):
    pivot = next(i for i,q in enumerate(source) if all(q))
    count = 0
    for perm in permutations(range(4)):
        ordered = [tuple(q[j] for j in perm) for q in source]
        for target_pivot in target:
            if not all(target_pivot):
                continue
            d = [t*pow(v,-1,prime)%prime for t,v in zip(target_pivot,ordered[pivot])]
            image = {normalize(tuple(z*y for z,y in zip(q,d)),prime) for q in ordered}
            count += image==set(target)
    return count


def invariants(u,v,a):
    z = v[0,0]*v[1,1]*u[0,1]*u[1,0]/(v[0,1]*v[1,0]*u[0,0]*u[1,1])
    left = u[0,1]*v[0,0]/(u[1,1]*v[0,1])
    right = u[0,0]*v[0,1]/(u[1,0]*v[0,0])
    return z,right*a[0,1]/a[0,0],left*a[1,0]/a[0,0],left*right*a[1,1]/a[0,0]


if __name__ == '__main__':
    u = sp.Matrix([[1,2],[3,5]])
    v = sp.Matrix([[2,3],[5,7]])
    a = sp.Matrix([[2,1],[5,3]])
    w = gram(u,v,a,other_completion(a,*cross_blocks(u,v)))
    pts = modular_points(w,P)
    def odd_points(matrix):
        def reduce(x):
            return int(x.p)*pow(int(x.q),-1,P)%P
        rel = [sp.Matrix(4,4,col).applyfunc(reduce) for col in matrix.columnspace()]
        odd = [kernel(sp.Matrix.vstack(*(sp.Matrix(1,4,q)*r for r in rel))) for q in pts]
        assert all(kernel(sp.Matrix.vstack(*((r*sp.Matrix(o)).T for r in rel)))==q
                   for q,o in zip(pts,odd))
        return odd
    configs = [framed(q) for q in [pts,odd_points(w),odd_points(w.T)]]
    assert [[number_of_maps(x,y) for y in configs] for x in configs] == [[4,0,0],[0,4,0],[0,0,4]]
    e = sp.symbols('e:4',nonzero=True)
    f = sp.symbols('f:4',nonzero=True)
    t = sp.symbols('t',nonzero=True)
    def thin(i,j):
        return e[i]*f[j]*t**(1 if (i^j)==2 else -1)
    us = sp.Matrix(2,2,lambda i,j:u[i,j]*thin(i,2+j))
    vs = sp.Matrix(2,2,lambda i,j:v[i,j]*thin(2+i,j))
    aa = sp.Matrix(2,2,lambda i,j:a[i,j]*e[i]*f[i]*e[j]*f[j^1])
    assert all(sp.cancel(i-j)==0 for i,j in zip(invariants(us,vs,aa),invariants(u,v,a)))
    print('Projective map counts diag(4,4,4); all marked rational invariants preserved.')
