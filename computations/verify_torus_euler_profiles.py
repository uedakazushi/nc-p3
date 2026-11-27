"""Exhaust all adjacent rank profiles for actions with two consecutive weights per slot."""
from itertools import combinations_with_replacement, permutations, product
from math import prod

K = range(4)
PERMS = [(p, (-1)**sum(p[i] > p[j] for i in K for j in range(i+1,4)))
         for p in permutations(K)]


def feasible(k, start, h):
    def dimension(i, weight):
        first, second = (4-k[i], k[i]), (4-k[(i+1)%4], k[(i+1)%4])
        return sum(first[u]*second[v] for u,v in product((0,1),repeat=2) if u+v == weight)
    for a in range(7):
        for c in range(7-a):
            row = (a, 6-a-c, c)
            if all(n <= min(dimension(start,j), dimension((start+2)%4,h-j))
                   for j,n in enumerate(row)):
                yield row


def opposite(row, h):
    return tuple(row[h-j] if 0 <= h-j <= 2 else 0 for j in range(3))


def determinant(k, h, relations, t, y):
    rows = []
    for i in K:
        prev = k[(i+3)%4]
        terms = (1+y*t**h, -(4-k[i]+k[i]*t),
                 sum(v*t**j for j,v in enumerate(relations[i])),
                 -t**(h-1)*((4-prev)*t+prev))
        rows.append([terms[(j-i)%4]*(y if j<i else 1) for j in K])
    return sum(sign*prod(rows[i][perm[i]] for i in K) for perm,sign in PERMS)


def profiles():
    # Keep the normalized scalar slots at weight zero, and retain a marked start up to recording rotations.
    result = set()
    for h in (1,2,3):
        factors = [v for v in combinations_with_replacement(range(5),4) if sum(v)==4*h]
        lookup = {prod(1-7*11**j for j in v):v for v in factors}
        assert len(lookup) == len(factors)
        for k in product(K,repeat=4):
            if k != min(k[i:]+k[:i] for i in K):
                continue
            for a,b in product(feasible(k,0,h), feasible(k,1,h)):
                r = (a,b,opposite(a,h),opposite(b,h))
                exponents = lookup.get(determinant(k,h,r,11,7))
                if exponents is None:
                    continue
                # Degree at most 4 in y and at most 16 in t on these charts.
                for t,y in product(range(1,18),range(5)):
                    assert determinant(k,h,r,t,y) == prod(1-y*t**j for j in exponents)
                assert sum(k) == 4*h
                result.add((k,h,a,b,exponents))
    return result


if __name__ == '__main__':
    assert profiles() == {
        ((0,1,0,3),1,(5,1,0),(5,1,0),(1,1,1,1)),
        ((0,1,2,1),1,(4,2,0),(2,4,0),(0,0,2,2)),
        ((1,1,1,1),1,(3,3,0),(3,3,0),(0,0,0,4)),
        ((0,3,2,3),2,(2,4,0),(0,4,2),(1,1,3,3)),
        ((1,1,3,3),2,(3,2,1),(0,6,0),(2,2,2,2)),
        ((1,1,3,3),2,(4,2,0),(1,4,1),(2,2,2,2)),
        ((2,2,2,2),2,(1,4,1),(1,4,1),(0,0,4,4)),
        ((3,3,3,3),3,(0,3,3),(0,3,3),(0,4,4,4)),
    }
    print('Eight determinant profiles before the consecutive-pure exclusion.')
