"""A finite free module test without an a priori all-degree Hilbert bound."""
from math import comb

import numpy as np
from flint import nmod_mat

from verify_klein_period_four_center import artinian_quotient, commuting_space, paths
from verify_parallel_invertible_rectangles import D, K, completed, matrix, product, quadric, s
from verify_parallel_third_rectangle import third_rectangle
from verify_quartic_euler_derivations import check_relations

PRIME = 1009


def residue(v):
    num,den = s.fraction(v)
    return int(num)*pow(int(den),-1,PRIME)%PRIME


def relations(w):
    return [[(h,[v[i] if i^j==h else s.S.Zero for i,j in product(K,repeat=2)])
             for h in K for v in matrix(w,p,h).columnspace()] for p in K]


def raw_dimensions(rel,start,limit):
    blocks = {(p,h):[[residue(row[4*i+(i^h)]) for i in K] for g,row in rs if g==h]
              for p,rs in enumerate(rel) for h in K}
    dims = [[1,0,0,0]]
    previous = {}
    for n in range(1,limit+1):
        current,arrows = [],{}
        for h in K:
            sizes = [dims[-1][h^j] for j in K]
            off = np.cumsum([0]+sizes)
            size = sum(sizes)
            rows = []
            if n>1:
                for b in K:
                    g = h^b
                    for col in blocks[(start+n-2)%4,b]:
                        block = np.zeros((dims[-2][g],size),dtype=np.int64)
                        for i,v in enumerate(col):
                            j = i^b
                            block[:,off[j]:off[j+1]]=v*previous[g,i]%PRIME
                        rows.extend(block.tolist())
            M = nmod_mat(rows,PRIME) if rows else nmod_mat(0,size,PRIME)
            reduced,rank = M.rref()
            dim = size-rank
            current.append(dim)
            projection = np.zeros((size,dim),dtype=np.int64)
            if rank:
                r = np.array(reduced.tolist()[:rank],dtype=np.int64)
                pivots = r.astype(bool).argmax(axis=1)
                free = np.delete(np.arange(size),pivots)
                projection[pivots,:]=(-r[:,free])%PRIME
            else:
                free = list(range(size))
            projection[free,:]=np.eye(dim,dtype=np.int64)
            for j in K:
                arrows[h^j,j]=projection[off[j]:off[j+1]]
        assert sum(current)==comb(n+3,3)
        previous = arrows
        dims.append(current)


def main():
    samples = []
    for shape in (76,112):
        t,b = s.Integer(2),s.Integer(3)
        c = (t**2-b)/(1-b)
        d = s.solve(quadric(shape,t,b,c,7,1,5,D),D)[0]
        samples.append(completed(shape,t,b,c,7,1,5,d))
    samples += [third_rectangle(1,2,2),third_rectangle(2,3,5)]
    for i,W in enumerate(samples):
        rel = relations(W)
        if i>=2:
            for p in K:
                raw_dimensions(rel,p,14)
        basis,img = paths(rel,6)
        assert [commuting_space(basis,img,h,0).cols for h in K]==[3,2,0,0]
        artinian_quotient(rel,basis,img)
        check_relations(rel,basis,img)
    print('PASS: all four quartic cuts and Euler derivation ranks, two independent residue border certificates')


if __name__=='__main__':
    main()
