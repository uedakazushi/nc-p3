"""One-sided missing colors: free paths, consecutive pure relations and cubic slice quotients."""
from functools import lru_cache
from itertools import product
from pathlib import Path
from subprocess import check_output, run
from tempfile import TemporaryDirectory

K = range(4)
SURVIVORS = (0x10330ffffff,0x660666ff9ff,0x66066ffffff,0x110120ffffff,0x110330ffffff,
             0x111111ffffff,0x330330ff0ff3,0x330330ff3ff3,0x330330ffffff,0x330cc0ffffff)


def decode(mask):
    return frozenset((i,j,j^e,i^e) for i,j,e in product(K,repeat=3) if i!=j and mask>>(12*e+4*((i^j)-1)+j)&1)


def encode(words):
    return sum(1<<(12*(i^l)+4*((i^j)-1)+j) for i,j,k,l in words)


def orbit(w):
    result = set()
    for words in (w,[k[2:]+k[:2] for k in w],[(b,a,d,c) for a,b,c,d in w],[(d,c,b,a) for a,b,c,d in w]):
        for pi in ((0,1,2,3),(0,1,3,2)):
            for s,h in product(K,(0,1)):
                result.add(encode((pi[a]^s,pi[b]^s,pi[c]^s^h,pi[d]^s^h) for a,b,c,d in words))
    return result


def rank(p,b):
    return (0 if b==0 else 2) if p%2==0 else (2 if b<2 else 1)


@lru_cache(None)
def lines(mask):
    rows = [(mask >> (4*i))&15 for i in K]
    active_columns = rows[0]|rows[1]|rows[2]|rows[3]
    parent = list(range(8))
    def root(i):
        while parent[i]!=i:
            i = parent[i]
        return i
    for i,j in product(K,repeat=2):
        if rows[i] and active_columns>>j&1 and not rows[i]>>j&1:
            parent[root(i)] = root(j+4)
    return rows,tuple(root(i) for i in K)


def upper(mask,r,indices):
    rows,classes = lines(mask)
    return min(r,len({classes[i] for i in indices if rows[i]}))


def admissible(w, block_rank=rank, row_bound=upper):
    mats = [[sum(1<<(4*t[p]+t[(p+2)%4]) for t in w if t[p]^t[(p+1)%4]==b) for b in K] for p in K]
    free = [[[(b^i) for b in K if not lines(mats[p][b])[0][i]] for i in K] for p in K]
    pure = [set() for _ in K]
    for p,b,i in product(K,repeat=3):
        r = block_rank(p,b)
        if r and row_bound(mats[p][b],r,[j for j in K if i!=j])<r:
            pure[p].add((i,i^b))
    if any(j==k for p in K for i,j in pure[p] for k,l in pure[(p+1)%4]):
        return False
    for p,i in product(K,repeat=2):
        x = free[p][i]
        z = [j for j in K if i in free[(p-1)%4][j]]
        lower = 0
        for g in K:
            d = len(x)*len(z)
            for c in z:
                b = c^g
                d -= row_bound(mats[(p+1)%4][b],block_rank(p+1,b),x)
            for a in x:
                b = a^g
                d -= row_bound(mats[(p+2)%4][b],block_rank(p+2,b),[c^b for c in z])
            lower += max(d,0)
        if lower>=2:
            return False
    return True


def main():
    source = Path(__file__).with_name('klein_single_missing_supports.cpp')
    with TemporaryDirectory() as temp:
        binary = Path(temp)/'supports'
        run(['c++','-std=c++17','-O3',str(source),'-o',str(binary)],check=True)
        masks = {int(m) for m in check_output([str(binary)],text=True).splitlines()}
    assert len(masks)==2664
    surviving = []
    count = 0
    while masks:
        mask = min(masks)
        w = decode(mask)
        group = orbit(w)
        assert group<=masks
        masks -= group
        count += 1
        if admissible(w):
            surviving.append(mask)
    assert count==74 and tuple(surviving)==SURVIVORS
    print('PASS: both one-sided rank profiles, 2664 nonexponential supports, 74 orbits, ten surviving tests')


if __name__=='__main__':
    main()
