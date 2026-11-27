"""Verify the generic unmarked stabilizer and the shift periods on the negative component."""
from functools import reduce
from itertools import permutations,product
from operator import xor

from flint import fmpz_mat
import sympy as s
from verify_klein_elliptic_completion import gram_twist


def main():
    K = range(4)
    words = [q for q in product(K,repeat=4) if reduce(xor,q)==0]
    degree = fmpz_mat([[int(q[p]==i) for p in K for i in K] for q in words])
    hermite,transform = degree.hnf(transform=True)
    smith = degree.snf()
    rank = degree.rank()
    assert [int(smith[i,i]) for i in range(rank)]==[1]*11+[2]*2
    assert transform*degree==hermite and abs(int(transform.det()))==1
    relations = [[(i,int(e)) for i,e in enumerate(row) if e] for row in transform.tolist()[rank:]]
    prime = 1009
    A = s.Matrix([[5,2],[3,7]])
    W = gram_twist(A,-A,s.Integer(3),s.S.NegativeOne)[0]
    def residue(x):
        return int(x.p)*pow(int(x.q),-1,prime)%prime
    coeff = {q:residue(W[4*q[0]+q[1],4*q[2]+q[3]]) for q in words}
    assert all(coeff.values())
    matches = []
    for shift in K:
        for permutation in permutations((1,2,3)):
            phi = (0,)+permutation
            for i,j,k in product(K,repeat=3):
                offsets = (i,j,k,i^j^k)
                ratios = []
                for q in words:
                    m = tuple(phi[l]^t for l,t in zip(q,offsets))
                    m = m[shift:]+m[:shift]
                    ratios.append(coeff[m]*pow(coeff[q],-1,prime)%prime)
                if all(s.prod(pow(ratios[i],e,prime) for i,e in row)%prime==1 for row in relations):
                    matches.append((shift,phi,offsets))
    assert matches==[(0,tuple(K),(0,0,0,0)),(2,tuple(K),(0,1,0,1))]
    print('PASS: the generic negative stabilizer has four elements and the minimum period is two')


if __name__=='__main__':
    main()
