"""Lower bounds for the distinct support counts of all four projective configurations."""
from itertools import product
import random

from verify_klein_elliptic_completion import gram_twist
from verify_rectangle_component_smoothing import K,examples,matrix,point_algebra,s


def main():
    A = s.Matrix([[5,2],[3,7]])
    W,_,_ = gram_twist(A,A,s.Integer(3),s.S.One)
    positive = {q:W[4*q[0]+q[1],4*q[2]+q[3]] for q in product(K,repeat=4)}
    for tensor,lower in zip(examples()+[positive],(18,19,18,20)):
        for p in K:
            rel = [(h,col) for h in K for col in matrix(tensor,p,h).columnspace()]
            rng = random.Random(772)
            mats,_ = point_algebra(rel,[rng.randrange(1,30) for _ in K])
            f = sum((v*2**j for j,v in enumerate(mats)),mats[0]*0).charpoly()
            assert (f//f.gcd(f.derivative())).degree()==lower
    print('PASS: the three horizontal representatives and positive pairing have at least 18,19,18,20 distinct points at each vertex')


if __name__=='__main__':
    main()
