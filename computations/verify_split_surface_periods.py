"""Half-shift characters on the two intrinsic deformation spaces of the Bell surface."""
from itertools import product

import sympy as sp
from flint import nmod_mat

from verify_split_bell_locus import GROUP
from verify_split_inverse_windows import BITS, bell, bell_coordinates
from verify_split_surface_lifts import WORDS, U, V, split_tensor

PRIME = 109
IMAG = 33


def matrix(rows):
    return nmod_mat(rows, PRIME)


def select(a,rows,cols):
    return matrix([[int(a[i,j]) for j in cols] for i in rows])


def pivots(a):
    reduced,rank = a.rref()
    return [next(j for j in range(a.ncols()) if reduced[i,j]) for i in range(rank)]


def kernel(a):
    k,dim = a.nullspace()
    return select(k,range(k.nrows()),range(dim))


def flatten(w,start):
    return nmod_mat(16,16,[w.get(q[-start:]+q[:-start] if start else q,0) for q in WORDS],PRIME)


def verify_completion(h, expected):
    generic = split_tensor(h)
    kappa = -(U**2-1)*(V**2-1)/(4*U*V)
    for q,value in generic.items():
        sector = tuple(i for i in range(4) if q[i]<2)
        factor = {(0,1):kappa,(2,3):1/kappa}.get(sector,1)
        assert sp.cancel(generic[q[2:]+q[:2]]-factor*value) == 0
    def residue(v):
        a = sp.Rational(sp.sympify(v).subs({U:2,V:3}))
        return int(a.p)*pow(int(a.q),-1,PRIME)%PRIME
    w = {q:residue(v) for q,v in generic.items()}
    bases = [(kernel(m.transpose()).transpose(),kernel(m)) for m in (flatten(w,0),flatten(w,1))]
    for n in (3,4):
        words = [q for q in WORDS if sum(j<2 for j in q) == n]
        cols = [[int(v) for i,(l,r) in enumerate(bases) for v in (l*flatten({q:1},i)*r).entries()]
                for q in words]
        k = kernel(matrix(cols).transpose())
        if n == 3:
            g = matrix([[w.get(q[:i]+(t,)+q[i+1:],0) if q[i]==a else 0
                         for i,a,t in product(range(4),range(2),range(2,4))] for q in words])
        else:
            g = nmod_mat(len(words),0,PRIME)
        assert (k.ncols(),g.rank()) == ((20,16) if n==3 else (3,0))
        both = matrix([[int(g[i,j]) for j in range(g.ncols())]+[int(k[i,j]) for j in range(k.ncols())]
                       for i in range(len(words))])
        frame = select(both,range(len(words)),pivots(both))
        rs = pivots(frame.transpose())
        lift = select(frame,rs,range(frame.ncols())).inv()
        free = range(g.rank(),frame.ncols())
        z = select(frame,range(len(words)),free)
        positions = {q:i for i,q in enumerate(words)}
        if n == 3:
            signed = matrix([[(-1)**sum(q)*int(z[i,j]) for j in range(z.ncols())]
                             for i,q in enumerate(words)])
            image = lift*select(signed,rs,range(z.ncols()))
            assert frame*image == signed
            action = select(image,free,range(z.ncols()))
            identity = matrix([[int(i==j) for j in range(4)] for i in range(4)])
            assert (action-identity).rank() == 2
        plus = []
        for flip,sign in product(range(2),repeat=2):
            transformed = []
            for q in words:
                qt = q[2:]+q[:2]
                scalar = sp.prod((-IMAG,-IMAG,IMAG,IMAG)[i] for i,c in enumerate(q) if c<2)
                if flip:
                    scalar *= sp.prod((IMAG,IMAG,-IMAG,-IMAG)[i] for i,c in enumerate(qt) if c<2)
                    qt = tuple(c^1 for c in qt)
                if sign:
                    scalar *= (-1)**sum(qt)
                transformed.append([int(scalar)*int(z[positions[qt],j])%PRIME for j in range(z.ncols())])
            image = lift*select(matrix(transformed),rs,range(z.ncols()))
            assert frame*image == matrix(transformed)
            action = select(image,free,range(z.ncols()))
            identity = matrix([[int(i==j) for j in range(z.ncols())] for i in range(z.ncols())])
            assert action*action == identity
            plus.append(z.ncols()-(action-identity).rank())
        assert plus == ([2]*4 if n == 3 else expected)


def verify():
    # This second sample distinguishes a one-shift even after exchanging both axes.
    spec = list(map(sp.Integer,[1,-28,4,7]))
    for tensor in (bell(spec),bell([1/x for x in spec])):
        rotated = bell_coordinates({q:tensor[q[1:]+q[:1]] for q in BITS})
        assert not any(all(rotated[p[j]]*e[j]*spec[0] == rotated[p[0]]*e[0]*spec[j]
                           for j in range(4)) for p,e in GROUP)
    x,z = sp.Matrix([[0,1],[1,0]]),sp.diag(1,-1)
    for h,expected in zip((sp.eye(2),x,z,x*z),([1,3,1,1],[1,1,1,3],[3,1,1,1],[1,1,3,1])):
        verify_completion(h,expected)
    print("All sixteen corrected half-shifts and both possible one-shift markings checked.")


if __name__ == '__main__':
    verify()
