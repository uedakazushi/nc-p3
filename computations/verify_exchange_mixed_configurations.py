"""Projective map bounds at all four starts of the first mixed-corner sample."""
from collections import Counter
from itertools import combinations, product

import sympy as s

from verify_klein_elliptic_completion import gram_twist
from verify_klein_full_rank_configurations import partner
from verify_klein_generic_isomorphisms import normalize
from verify_klein_two_doublets import modular_points


def incidence_signatures(points, prime):
    planes = set()
    triples = []
    for I in combinations(range(len(points)),3):
        a,b,c = (points[i] for i in I)
        coeff = []
        for i in range(4):
            j,k,l = [n for n in range(4) if n!=i]
            det = (a[j]*(b[k]*c[l]-b[l]*c[k])-a[k]*(b[j]*c[l]-b[l]*c[j])
                   +a[l]*(b[j]*c[k]-b[k]*c[j]))
            coeff.append((-1)**i*det%prime)
        if not any(coeff):
            triples.append(I)
            continue
        inds = tuple(i for i,q in enumerate(points) if sum(x*y for x,y in zip(coeff,q))%prime==0)
        if len(inds)>3:
            planes.add(inds)
    histogram = Counter(map(len,planes))
    signatures = [tuple(sum(i in pl for pl in planes if len(pl)==n) for n in sorted(histogram))
                  + (sum(i in I for I in triples),) for i in range(len(points))]
    return histogram,signatures


def incidences(points, prime):
    histogram,signatures = incidence_signatures(points,prime)
    assert histogram=={4:52} and all(sig[-1]==0 for sig in signatures)
    return [sig[0] for sig in signatures]


def maps(source,target,prime):
    src,a = source
    dst,b = target
    if Counter(a)!=Counter(b):
        return 0
    def matrix(pts,I):
        return s.Matrix([pts[i] for i in I]).T
    def coordinates(inv,q):
        return tuple(int(v)%prime for v in inv*s.Matrix(q))
    for I in sorted(combinations(range(20),4),key=lambda I:s.prod(a.count(a[i]) for i in I)):
        A = matrix(src,I)
        if int(A.det())%prime:
            inv = A.inv_mod(prime)
            frame = [normalize(coordinates(inv,q),prime) for q in src]
            pivot = next((i for i,q in enumerate(frame) if all(q)),None)
            if pivot is not None:
                break
    else:
        raise AssertionError("source has no projective frame")
    count = 0
    for images in product(*[[j for j in range(20) if b[j]==a[i]] for i in I]):
        if len(set(images))!=4:
            continue
        A = matrix(dst,images)
        if not int(A.det())%prime:
            continue
        inv = A.inv_mod(prime)
        targets = [normalize(coordinates(inv,q),prime) for q in dst]
        for j,q in enumerate(targets):
            if b[j]!=a[pivot] or not all(q):
                continue
            diag = [y*pow(v,-1,prime)%prime for y,v in zip(q,frame[pivot])]
            count += {normalize(tuple(v*z for v,z in zip(pt,diag)),prime) for pt in frame}==set(targets)
    return count


if __name__ == "__main__":
    p = 1753
    A = s.Matrix([[2,3],[3,7]])
    B = s.Matrix([[3,5],[s.Rational(2,5),3]])
    w,_,n = gram_twist(A,B,s.S.One,s.Integer(2))
    configs = [modular_points(w,p)]
    for mat,transpose in [(w,False),(n,True),(w.T,False)]:
        configs.append(partner(configs[-1],mat,transpose,p))
    assert partner(configs[-1],n.T,True,p)==configs[0]
    flip = [4*j+i for i,j in product(range(4),repeat=2)]
    for pts,mat in zip(configs[1:],(n[flip,:],w.T,n.T[flip,:])):
        assert set(modular_points(mat,p))==set(pts)
    tagged = [(pts,incidences(pts,p)) for pts in configs]
    assert [[maps(x,y,p) for y in tagged] for x in tagged]==[
        [4,0,0,4],[0,4,4,0],[0,4,4,0],[4,0,0,4]]
    print("PASS: four reduced configurations and all projective map counts modulo 1753")
